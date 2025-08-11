import requests
from bs4 import BeautifulSoup
import pandas as pd
from time import sleep

# STEP 1: Ensembl API를 사용해서 REF/ALT allele 가져오기 (멀티 ALT 처리 포함)
def get_ref_alt_from_ensembl(rsid):
    url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
    response = requests.get(url)
    if response.status_code != 200:
        redirect_url = response.headers.get('location')
        if redirect_url and 'rs' in redirect_url:
            merged_rsid = redirect_url.split('/')[-1]
            print(f"🔁 [{rsid}] → 발견 ID: {merged_rsid}")
            return get_ref_alt_from_ensembl(merged_rsid)
        raise Exception(f"[{rsid}] Ensembl API 요청 실패: {response.status_code}")

    data = response.json()
    if 'merged' in data and data['merged']:
        merged_into = data['merged'][0]['id']
        print(f"🔁 [{rsid}] → 발견 ID: {merged_into}")
        return get_ref_alt_from_ensembl(merged_into)

    for mapping in data.get('mappings', []):
        if mapping['assembly_name'] in ['GRCh38', 'GRCh37']:
            alleles = mapping['allele_string'].split('/')
            if len(alleles) >= 2:
                ref = alleles[0].upper()
                alts = [a.upper() for a in alleles[1:]]
                return ref, alts, rsid

    raise Exception(f"[{rsid}] GRCh38 또는 GRCh37 기반의 REF/ALT 정보를 찾을 수 없음")

# STEP 2.1: 발견 ID 추적 (dbSNP 기준)
def resolve_merged_rsid_from_dbsnp(rsid):
    url = f"https://www.ncbi.nlm.nih.gov/snp/{rsid}"
    headers = {'User-Agent': 'Mozilla/5.0'}
    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        return rsid

    soup = BeautifulSoup(response.text, 'html.parser')
    link = soup.find('a', href=True, string=lambda s: s and s.lower().startswith('rs'))
    if link and 'merged into' in link.find_parent('div').text.lower():
        merged_id = link.text.strip()
        print(f"🔁 [dbSNP] {rsid} → 발견 ID: {merged_id}")
        return merged_id
    return rsid

# STEP 2.2: dbSNP ALT allele 모두 저장 + Ensembl 포함 여부 표시
def get_frequency_from_dbsnp(rsid, ref, alt_list, original_rsid=None):
    url = f"https://www.ncbi.nlm.nih.gov/snp/{rsid}"
    headers = {'User-Agent': 'Mozilla/5.0'}
    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        raise Exception(f"[{rsid}] dbSNP 페이지 요청 실패: {response.status_code}")

    soup = BeautifulSoup(response.text, 'html.parser')
    tables = soup.find_all('table')

    results = []
    alt_str = ",".join(alt_list)

    for table in tables:
        rows = table.find_all('tr')
        if not rows:
            continue

        header_cells = rows[0].find_all(['th', 'td'])
        headers = [cell.text.strip().upper() for cell in header_cells]
        print(f"[{rsid}] 헤더: {headers}")

        if 'REF ALLELE' in headers and 'ALT ALLELE' in headers:
            ref_index = headers.index('REF ALLELE')
            alt_index = headers.index('ALT ALLELE')
            pop_index = headers.index('POPULATION') if 'POPULATION' in headers else 0
            count_index = headers.index('SAMPLE SIZE') if 'SAMPLE SIZE' in headers else -1

            for row in rows[1:]:
                cols = [td.text.strip() for td in row.find_all('td')]
                if len(cols) <= max(ref_index, alt_index):
                    continue

                population = cols[pop_index] if pop_index >= 0 else "UNKNOWN"
                allele_number = cols[count_index] if count_index >= 0 else ""

                ref_parts = cols[ref_index].split('=')
                if len(ref_parts) == 2:
                    ref_allele = ref_parts[0].strip().upper()
                    try:
                        ref_freq = float(ref_parts[1].replace(',', '').strip())
                    except:
                        ref_freq = None
                else:
                    ref_allele = ""
                    ref_freq = None

                if ref_allele:
                    results.append({
                        'rsID': original_rsid or rsid,
                        'merged_from': rsid if original_rsid and original_rsid != rsid else '',
                        'population': population,
                        'allele': ref_allele,
                        'allele_type': "REF",
                        'frequency': ref_freq,
                        'allele_number': allele_number,
                        'alt_alleles': alt_str,
                        'is_in_ensembl': True
                    })

                alt_allele_entries = cols[alt_index].split(',')
                for entry in alt_allele_entries:
                    alt_parts = entry.split('=')
                    if len(alt_parts) == 2:
                        alt_allele = alt_parts[0].strip().upper()
                        try:
                            alt_freq = float(alt_parts[1].replace(',', '').strip())
                        except:
                            alt_freq = None
                    else:
                        continue

                    is_in_ensembl = alt_allele in alt_list
                    results.append({
                        'rsID': original_rsid or rsid,
                        'merged_from': rsid if original_rsid and original_rsid != rsid else '',
                        'population': population,
                        'allele': alt_allele,
                        'allele_type': "ALT",
                        'frequency': alt_freq,
                        'allele_number': allele_number,
                        'alt_alleles': alt_str,
                        'is_in_ensembl': is_in_ensembl
                    })
    return results

# STEP 3: rsID 보드 처리
def process_rsids(rsid_list):
    all_data = []
    errors = []

    for original_rsid in rsid_list:
        print(f"🔍 처리 중: {original_rsid}")
        try:
            merged_rsid_dbsnp = resolve_merged_rsid_from_dbsnp(original_rsid)
            ref, alts, dbsnp_rsid = get_ref_alt_from_ensembl(merged_rsid_dbsnp)
            result = get_frequency_from_dbsnp(dbsnp_rsid, ref, alts, original_rsid=original_rsid)
            for row in result:
                row['rsID'] = original_rsid
                row['merged_from'] = dbsnp_rsid if original_rsid != dbsnp_rsid else ''
            all_data.extend(result)
        except Exception as e:
            print(f"⚠️ 오류 발생: {e}")
            errors.append({'rsID': original_rsid, 'error': str(e)})
        sleep(1)

    return all_data, errors

# STEP 4: population summary (long format)
def summarize_by_population(data):
    df = pd.DataFrame(data)
    df['frequency'] = pd.to_numeric(df['frequency'], errors='coerce')
    df['allele_number'] = pd.to_numeric(df['allele_number'], errors='coerce')

    summary_rows = []
    grouped = df.groupby(['rsID', 'population'])

    for (rsid, population), group in grouped:
        row = {'rsID': rsid, 'population': population, 'sample_size': None}
        ref_row = group[group['allele_type'] == 'REF']
        alt_row = group[group['allele_type'] == 'ALT']

        if not ref_row.empty:
            row['REF_allele'] = ref_row.iloc[0]['allele']
            row['REF_freq'] = ref_row.iloc[0]['frequency']
            row['sample_size'] = ref_row.iloc[0]['allele_number']
            row['alt_alleles'] = ref_row.iloc[0].get('alt_alleles', '')

        if not alt_row.empty:
            row['ALT_allele'] = alt_row.iloc[0]['allele']
            row['ALT_freq'] = alt_row.iloc[0]['frequency']
            if row['sample_size'] is None:
                row['sample_size'] = alt_row.iloc[0]['allele_number']
            row['alt_alleles'] = alt_row.iloc[0].get('alt_alleles', '')

        if 'REF_freq' in row or 'ALT_freq' in row:
            summary_rows.append(row)

    return pd.DataFrame(summary_rows)

# STEP 5: 넓은 형식 요약 + ALT allele 이름 포함
def reshape_population_summary_wide(data):
    df = pd.DataFrame(data)
    df['frequency'] = pd.to_numeric(df['frequency'], errors='coerce')

    ref_df = df[df['allele_type'] == 'REF'][['rsID', 'population', 'allele', 'frequency']]
    alt_df = df[df['allele_type'] == 'ALT'][['rsID', 'population', 'allele', 'frequency']]

    ref_df = ref_df.rename(columns={'allele': 'REF_allele', 'frequency': 'REF_freq'})
    alt_df = alt_df.rename(columns={'allele': 'ALT_allele', 'frequency': 'ALT_freq'})

    merged = pd.merge(ref_df, alt_df, on=['rsID', 'population'], how='outer')

    ref_pivot = merged.pivot_table(index='rsID', columns='population', values='REF_freq', aggfunc='first')
    alt_freq_pivot = merged.pivot_table(index='rsID', columns='population', values='ALT_freq', aggfunc='first')
    alt_allele_pivot = merged.pivot_table(index='rsID', columns='population', values='ALT_allele', aggfunc='first')

    ref_pivot.columns = [f"{pop}_REF" for pop in ref_pivot.columns]
    alt_freq_pivot.columns = [f"{pop}_ALT" for pop in alt_freq_pivot.columns]
    alt_allele_pivot.columns = [f"{pop}_ALT_allele" for pop in alt_allele_pivot.columns]

    allele_info = merged.groupby('rsID').first().reset_index()[['rsID', 'REF_allele', 'ALT_allele']]

    final_df = allele_info.set_index('rsID').join([
        ref_pivot, alt_freq_pivot, alt_allele_pivot
    ]).reset_index()

    return final_df

# STEP 6: 저장 함수
def save_outputs(data, errors,
                 excel_file="combined_frequencies.xlsx",
                 csv_file="combined_frequencies.csv"):

    df = pd.DataFrame(data)
    df_summary = summarize_by_population(data)
    df_wide = reshape_population_summary_wide(data)

    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name="Allele_Frequencies", index=False)
        df_summary.to_excel(writer, sheet_name="Population_Summary", index=False)
        df_wide.to_excel(writer, sheet_name="Population_Wide", index=False)
        if errors:
            df_errors = pd.DataFrame(errors)
            df_errors.to_excel(writer, sheet_name="Errors", index=False)

    df.to_csv(csv_file, index=False)
    df_summary.to_csv("population_summary.csv", index=False)
    df_wide.to_csv("population_summary_wide.csv", index=False)
    if errors:
        df_errors.to_csv("frequency_errors.csv", index=False)

    print("✅ 저장 완료:")
    print(f"  - {excel_file}")
    print(f"  - {csv_file}")
    print(f"  - population_summary.csv")
    print(f"  - population_summary_wide.csv")

# STEP 7: 실행 예시
if __name__ == "__main__":

    rsid_list = [

        "rs10073425", "rs10152544", "rs10174949", "rs10214273", "rs1059513", "rs10738626", "rs10790275", "rs10791824",
        "rs10796303", "rs10822037", "rs10833", "rs10836538", "rs10888499", "rs10988863", "rs10995245", "rs10995251",
        "rs10995255", "rs10995256", "rs11033603", "rs11171739", "rs11205006", "rs112111458", "rs11216206", "rs11236791",
        "rs11236809", "rs11236813", "rs112385344", "rs114503346", "rs115161931", "rs117137535", "rs11786685",
        "rs11811788", "rs118162691", "rs11857092", "rs12123821", "rs12126142", "rs12133641", "rs12138773", "rs12144049",
        "rs1214598", "rs12153855", "rs12244238", "rs12251307", "rs12295535", "rs12565349", "rs12586305", "rs12634229",
        "rs12743520", "rs1295686", "rs13015714", "rs13097010", "rs13403179", "rs1358175", "rs1409123", "rs1438673",
        "rs1444418", "rs1444789", "rs146527530", "rs148161264", "rs149199808", "rs150979174", "rs1665050", "rs16862519",
        "rs16999165", "rs17132590", "rs17368814", "rs17371133", "rs17389644", "rs176095", "rs17881320", "rs183884396",
        "rs1861246", "rs187080438", "rs187325802", "rs188069315", "rs188720898", "rs189163698", "rs2050190",
        "rs2075943", "rs2164983", "rs2212434", "rs2227491", "rs2259735", "rs2271404", "rs2272128", "rs2415269",
        "rs2426500", "rs2542147", "rs2766664", "rs280024", "rs28383323", "rs28383330", "rs28406364", "rs28520436",
        "rs28558565", "rs2897442", "rs2967677", "rs301804", "rs3091307", "rs3099143", "rs3125788", "rs3126085",
        "rs3208007", "rs34215892", "rs34290285", "rs35073649", "rs35570272", "rs35766269", "rs3757723", "rs3848669",
        "rs3853601", "rs3862469", "rs3864302", "rs3947727", "rs41268896", "rs41293876", "rs4131280", "rs4247364",
        "rs4262739", "rs4312054", "rs4532376", "rs45599938", "rs45605540", "rs4574025", "rs4705908", "rs4706020",
        "rs471144", "rs4713555", "rs4722404", "rs4759228", "rs4796793", "rs479844", "rs4821544", "rs4821569",
        "rs4845373", "rs4906263", "rs5005507", "rs538763482", "rs56101042", "rs56302621", "rs56308324", "rs5743614",
        "rs59039403", "rs593982", "rs6023002", "rs6062486", "rs61776548", "rs61815704", "rs61816766", "rs61865882",
        "rs61878692", "rs62193132", "rs629326", "rs6461503", "rs659529", "rs6661961", "rs6720763", "rs675531",
        "rs67766926", "rs6785012", "rs6808249", "rs6996614", "rs7000782", "rs705699", "rs7110818", "rs7127307",
        "rs7130588", "rs7147439", "rs71625130", "rs72702813", "rs72702900", "rs72823628", "rs72925996", "rs72943976",
        "rs73018933", "rs75024669", "rs7542147", "rs759382", "rs7613051", "rs7701967", "rs7717955", "rs7815944",
        "rs7843258", "rs7857407", "rs78914480", "rs7927894", "rs7936323", "rs79497729", "rs80199341", "rs8086340",
        "rs821429", "rs847", "rs848", "rs859723", "rs878860", "rs891058", "rs909341", "rs9275218", "rs9368677",
        "rs9469099", "rs952558", "rs9540294", "rs9540298", "rs9864845", "rs989437", "rs9911533", "rs9923856",
        "rs1003878", "rs1024161", "rs1033500", "rs1063355", "rs10760706", "rs1077393", "rs10807113", "rs10876864",
        "rs10947262", "rs1107345", "rs11155700", "rs11752643", "rs11759611", "rs12177980", "rs12183587", "rs12202737",
        "rs12213837", "rs1270942", "rs13199787", "rs13729", "rs1413901", "rs142986308", "rs16898264", "rs1701704",
        "rs17429444", "rs17500468", "rs1794282", "rs1860545", "rs1980493", "rs2009345", "rs2010259", "rs2051549",
        "rs2069408", "rs2070600", "rs2071800", "rs2072633", "rs2076530", "rs2076537", "rs2137497", "rs2155219",
        "rs2187668", "rs2216164", "rs2239804", "rs2269426", "rs2292239", "rs2301271", "rs231726", "rs231735",
        "rs231775", "rs231804", "rs2395162", "rs2395163", "rs2395174", "rs2395175", "rs2395182", "rs2442749",
        "rs2476601", "rs2647012", "rs2647050", "rs2856717", "rs2856718", "rs2856725", "rs2858305", "rs2858331",
        "rs2858332", "rs2859078", "rs304650", "rs3096851", "rs3096866", "rs3104404", "rs3104405", "rs3115553",
        "rs3115573", "rs3116504", "rs3117099", "rs3118470", "rs3129871", "rs3129890", "rs3129939", "rs3129943",
        "rs3129963", "rs3130315", "rs3130320", "rs3130340", "rs3135353", "rs3763309", "rs3763312", "rs377763",
        "rs3789129", "rs3817973", "rs389883", "rs389884", "rs405875", "rs4147359", "rs4151657", "rs437179", "rs4424066",
        "rs470138", "rs494620", "rs547077", "rs547261", "rs574087", "rs6457536", "rs6457617", "rs652888", "rs653178",
        "rs6901084", "rs6903130", "rs6910071", "rs6911628", "rs6935051", "rs6935269", "rs6941112", "rs694739",
        "rs705708", "rs706779", "rs707928", "rs7192", "rs7453920", "rs7500151", "rs7682241", "rs7682481", "rs773107",
        "rs7745656", "rs7756516", "rs7758736", "rs7775397", "rs805294", "rs805303", "rs8111", "rs926169", "rs9267947",
        "rs9268132", "rs9268368", "rs9268384", "rs9268528", "rs9268530", "rs9268542", "rs9268615", "rs9268832",
        "rs9275224", "rs9275524", "rs9275572", "rs9275659", "rs9275686", "rs9275698", "rs9276435", "rs9357152",
        "rs9368713", "rs9397624", "rs9405090", "rs9461799", "rs9479482", "rs972099", "rs9864529", "rs10036748",
        "rs10046127", "rs1008953", "rs1042636", "rs1043483", "rs10484554", "rs1056198", "rs10737548", "rs1076160",
        "rs1077492", "rs10782001", "rs10789285", "rs10794648", "rs10832027", "rs10865331", "rs10888501", "rs10889668",
        "rs10960680", "rs11059675", "rs11065979", "rs11121129", "rs11209026", "rs11465802", "rs115324207", "rs11593576",
        "rs116432905", "rs11687879", "rs11746443", "rs11757367", "rs11795343", "rs11961853", "rs12042824", "rs12188300",
        "rs12191877", "rs12199223", "rs12206377", "rs1250544", "rs1250546", "rs1250566", "rs12524487", "rs12564022",
        "rs12580100", "rs12586317", "rs1265181", "rs12720356", "rs12790634", "rs12884468", "rs13437088", "rs1400473",
        "rs1473247", "rs147965700", "rs149912748", "rs1576", "rs1581803", "rs16895575", "rs16899661", "rs16949",
        "rs171329", "rs17177618", "rs17185076", "rs17259252", "rs17272796", "rs17728338", "rs181359", "rs1892497",
        "rs194675", "rs1967", "rs1975974", "rs1990760", "rs20541", "rs2066807", "rs2066808", "rs2066818", "rs2082412",
        "rs210192", "rs2111485", "rs2145623", "rs2179920", "rs2201841", "rs2229092", "rs2230653", "rs2233278",
        "rs2240804", "rs2256594", "rs2276405", "rs2328530", "rs2395029", "rs240993", "rs2451258", "rs2517600",
        "rs2523710", "rs2546890", "rs2675662", "rs2688608", "rs2700979", "rs2700984", "rs27524", "rs2778031",
        "rs280497", "rs28366363", "rs2836747", "rs2844579", "rs2844651", "rs28512356", "rs2853694", "rs2857597",
        "rs2857602", "rs28998802", "rs29261", "rs2944542", "rs3094165", "rs3116807", "rs3129207", "rs3129269",
        "rs3129817", "rs3129882", "rs3130192", "rs3130455", "rs3132572", "rs3174808", "rs3184504", "rs3213094",
        "rs33980500", "rs34115245", "rs34172843", "rs34394770", "rs35960711", "rs367254", "rs3729508", "rs3730013",
        "rs3747517", "rs3778620", "rs3782886", "rs3794765", "rs3802826", "rs397081", "rs4074995", "rs408036",
        "rs4085613", "rs4112788", "rs41268474", "rs4148871", "rs416603", "rs4242369", "rs4245080", "rs4313034",
        "rs4406273", "rs453779", "rs458017", "rs4649203", "rs465969", "rs4664464", "rs4713466", "rs4713614",
        "rs4795067", "rs4845454", "rs4908343", "rs4921493", "rs492602", "rs495337", "rs5063", "rs56095701",
        "rs57550632", "rs582757", "rs610037", "rs610604", "rs61839660", "rs61907765", "rs643177", "rs6444895",
        "rs6457109", "rs6457702", "rs6545930", "rs6590334", "rs6596086", "rs6677595", "rs671", "rs675640", "rs6759003",
        "rs6760912", "rs6809854", "rs6860806", "rs6870828", "rs6894567", "rs6903989", "rs694764", "rs696", "rs702873",
        "rs72866766", "rs72896150", "rs73127695", "rs736801", "rs740884", "rs744487", "rs75105906", "rs7548511",
        "rs76337351", "rs7637230", "rs76462670", "rs7665090", "rs7709212", "rs7720046", "rs7745603", "rs7769061",
        "rs78233367", "rs7993214", "rs8005252", "rs8016947", "rs8128234", "rs842625", "rs842636", "rs86567", "rs892085",
        "rs9260734", "rs9260740", "rs9266630", "rs9295676", "rs9295895", "rs9304742", "rs9348718", "rs9348841",
        "rs9366724", "rs9368611", "rs9380326", "rs9393967", "rs9393991", "rs9394026", "rs9400467", "rs9468487",
        "rs9504361", "rs9533962", "rs984971", "rs10155912", "rs10200159", "rs10250629", "rs1031034", "rs1042602",
        "rs1043101", "rs10768122", "rs10774624", "rs10986311", "rs11021232", "rs11079035", "rs11203203", "rs1126809",
        "rs1129038", "rs11966200", "rs12203592", "rs12421615", "rs12482904", "rs12771452", "rs12973771", "rs13076312",
        "rs13136820", "rs13208776", "rs13227879", "rs1393350", "rs1417210", "rs1464510", "rs148136154", "rs1635168",
        "rs16843742", "rs16872571", "rs17128310", "rs1805007.2", "rs2017445", "rs210124", "rs2122476", "rs2236313",
        "rs2247314", "rs2274195", "rs229527", "rs2304206", "rs231725", "rs2456973", "rs251464", "rs2687812",
        "rs28366353", "rs301807", "rs3127197", "rs34346645", "rs35095897", "rs35161626", "rs35860234", "rs3757247",
        "rs3806156", "rs3814231", "rs3823355", "rs41342147", "rs4268748", "rs4308124", "rs4409785", "rs4766578",
        "rs4807000", "rs4822024", "rs4908760", "rs561079", "rs56207241", "rs59374417", "rs6012953", "rs6059655",
        "rs638893", "rs6510827", "rs6583331", "rs6679677", "rs6902119", "rs71508903", "rs7188793", "rs72928038",
        "rs7758128", "rs78521699", "rs8083511", "rs8192917", "rs853308", "rs870355", "rs9271597", "rs9380143",
        "rs9468925", "rs9611565", "rs968567", "rs982204", "rs9851967", "rs9926296"



    ]
    rsid_list2 = ["rs10833", "rs982204", "rs9851967", "rs9926296"]
    results, errors = process_rsids(rsid_list2)
    save_outputs(results, errors)
