# module imports
import requests
import json
import time
import pandas as pd
import urllib.parse
import urllib.request
from xml.etree import ElementTree
from datetime import datetime
import generate_database

# uncomment the below line to provide a hard-coded json file
# json_file = 'pmids_20.json'
# get the current date-time to create output file name
date = datetime.now().strftime("%d%d%Y_%H_%M_%S")


def load_json(file_name):
    """
    This function loads the input json file and returns the
    json data
    :param file_name: the JSON file that needs to be loaded
    :return: the JSON data
    """
    json_file_content = {}
    with open(file_name) as json_obj:
        json_file_content = json_obj.read().split('\n')
    return json_file_content


def parse_json(json_file_content):
    """
    This function extracts the PMIDs from the json data.
    :param json_file_content: the json data from which PMIDs need to be extracted
    :return: the list of PMIDs
    """
    pmids = []
    for json_data in json_file_content:
        if len(json_data) != 0:
            pmids.append(json.loads(json_data)['docId'])
    return pmids


def extract_pmids():
    """
    This method extracts the PMIDs using PubMed API
    :return: list of PMIDs
    """
    # the url to get PMIDs from PubMed API
    pubmed_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=Abiotic+AND+Stimulus+AND+Plant[Mesh]"

    # the payload for the HTTP POST request
    pubmed_payload = {}
    headers = {
        'Cookie': 'ncbi_sid=F9FA0057F1DE9983_815276SID'
    }

    # send HTTP post request and extract the json formatted output
    response = requests.request("POST", pubmed_url, headers=headers, data=json.dumps(pubmed_payload), verify=False)

    pubmed_result = response.text.encode('utf8')

    # extract only the PMIDs from the JSON response
    tree = ElementTree.fromstring(pubmed_result)
    ids = tree.findall("./IdList/Id")
    pmids = [pmid.text for pmid in ids]
    print('PMIDs: ', pmids)
    # return the result
    return pmids


def extract_pmids_from_json(file_name):
    """
    This method extracts the PMIDs from the JSON file
    :return: list of PMIDs
    """
    # invoke the function that parses the json file
    # and returns the list of pmids
    pmids = parse_json(file_name)
    return pmids


def extract_pgenn_response(pmids):
    """
    This method extracts the PGenn API for each PMID
    :param pmids: the list of PMIDs extracted using TextPresso API
    :return: the list containing pgenn response for each PMID
    """
    # initialize empty list
    pgenn_response = []
    # for this API, payload will be empty
    payload = {}
    headers = {
        'Content-Type': 'application/json'
    }
    # iterate through each pmid in pmids list and generate url for PGenn API
    # send HTTP request with this url
    for pmid in pmids:
        print("Extracting PGenn response for PMID: {0}...".format(pmid))
        pgenn_url = "https://research.bioinformatics.udel.edu/itextmine/api/search/pmid/pgenn/medline.aligned/" + pmid + "/"
        response = requests.request("GET", pgenn_url, headers=headers, data=payload, verify=False)
        # append the JSON formatted response for the url in the list
        pgenn_response.append(response.json())
        # give a delay between each HTTP requests so that multiple HTTP requests do not crash the server
        time.sleep(1) # increased from 0.5 -> 1 due to HTTP connection exception
        print("Done")
    # finally, return the list containing PGenn responses for each PMID
    return pgenn_response


def get_gene_uniprot(pgenn_response):
    """
    This method extracts the gene name and UniprotID for each PMID and generates
    a dictionary containing PMID as key and a tuple of gene name and UniprotID as values.
    If there are multiple genes the value for a PMID key will contain a list of tuples
    :param pgenn_response: the list containing pgenn response
    :return: a dictionary with key as PMID and value as list of tuples of the form (gene_name, uniprot_id).
    """
    # initialize empty dictionary
    pgenn_results = {}

    url = 'https://www.uniprot.org/uploadlists/'

    # iterate through each response from the pgenn response list
    for response in pgenn_response:
        # for certain PMIDs, there will no response, hence it will return an empty list
        if len(response) != 0:
            # if the response is not empty, extract the first element since all the responses
            # will be a list containing single element
            response = response[0]
            # iterate through each 'entity' attribute in the response
            for entity in response['entity']:
                # the 'docId' attribute holds the PMID for that response
                # 'name' attribute holds gene name and 'uniprot' attribute holds the UniprotID
                # in each 'entity' attributes in the JSON response
                # thus, for each 'entity' attribute, create a list of tuples and add it to the
                # PMID key in dictionary
                if response['docId'] not in pgenn_results.keys():
                    # if the PMID is not already present in dictionary, initialize it with an empty list
                    pgenn_results[response['docId']] = []

                # now, get the correct entrez ID for the given uniprot ID
                params = {
                    'from': 'ACC',
                    'to': 'P_ENTREZGENEID',
                    'format': 'tab',
                    'query': entity['uniprot']
                }
                print('Uniprot ID: ', entity['uniprot'])
                data = urllib.parse.urlencode(params)
                data = data.encode('utf-8')
                req = urllib.request.Request(url, data)
                with urllib.request.urlopen(req) as f:
                    result = f.read()
                result = result.decode('utf8').split('\n')[1].split('\t')
                if len(result) == 2:
                    entrez = result[1]
                else:
                    entrez = ""
                print('Retrieved Entrez ID: ', entrez)
                # now, append to the list for that specific PMID, a tuple of name and UniprotID
                pgenn_results[response['docId']].append((entity['name'], entity['uniprot'], entrez))
    # finally, return the dictionary
    return pgenn_results


def extract_pubtator_response(pmids):
    """
    This function extracts the response from Pubtator API for each PMID
    :param pmids: the list of PMIDs
    :return: the list containing HTTP responses from Pubtator API for each PMID
    """
    # initialize empty list to hold the HTTP responses
    pubtator_responses = []
    # for this API, payload is empty
    payload = {}
    headers = {
        'Cookie': 'ncbi_sid=F9FA0057F1DE9983_815276SID'
    }
    # iterate through each PMID, and generate the corresponding pubtator url for the PMID
    for pmid in pmids:
        print("Extracting Pubtator response for PMID: {0}...".format(pmid))
        pubtator_url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmids=" + pmid + "&concepts=gene"
        # send HTTP GET request with this url and append to the list, the JSON formatted output
        response = requests.request("GET", pubtator_url, headers=headers, data=payload, verify=False)
        pubtator_responses.append(response.json())
        # give a delay of 0.5s to avoid the server getting crashed due to multiple HTTP requests
        time.sleep(1)
        print("Done.")
    # finally, return the list containing the HTTP responses
    return pubtator_responses


def get_gene_entrez_id(pubtator_responses):
    """
    This function parses the HTTP responses to extract only the gene name and EntrezID
    for each PMID. In this scenario as well, we will have multiple genes and IDs. Hence the results
    will contain a list of tuples of the form (gene_name, entrez_id) for each PMID
    :param pubtator_responses: the list of HTTP responses from Pubtator API
    :return: the dictionary containing the parsed response
    """
    # initialize an empty dictionary to hold the parsed result
    pubtator_results = {}
    # iterate through each response in the response list
    for response in pubtator_responses:
        # the Gene Name and EntrezID is present in a deep nested child node
        # for each response, extract the 'passages' node. From this, extract
        # the 'annotations' node. This contains 'infons' node. This node contains 'type'
        # attribute which specifies whether its a Gene or not. If its a gene, then get the Gene
        # name and EntrezID. The 'id' attribute in response will hold the PMID. The 'annotations' node
        # contains 'text' which has the Gene Name and the 'identifier' attribute of 'infons' node contains the
        # EntrezID
        for annot in response['passages']:
            # iterate through passages, and then through each annotations
            for ann in annot['annotations']:
                # if the type attribute of infons node in each annotation is a Gene,
                if ann['infons']['type'] == 'Gene':
                    # if the 'id' attribute (PMID) is not already present in the dictionary,
                    if response['id'] not in pubtator_results.keys():
                        # initialize an empty dictionary for the specific id (PMID)
                        pubtator_results[response['id']] = []
                    # now, append to the list, a tuple of Gene Name and EntrezID corresponding to the PMID (key)
                    pubtator_results[response['id']].append((ann['text'], ann['infons']['identifier']))
        # finally, return the dictionary holding the parsed responses
        return pubtator_results


def get_annotations(pmids):
    """
    This method extracts annotations for each PMID
    :param pmids: the list of PMIDs
    :return: the list containing HTTP responses with annotation information for each PMID
    """
    # initialize an empty list to store the HTTP responses
    annotations = []
    # in this API, payload is empty
    payload = {}
    headers = {
        'Cookie': 'JSESSIONID=EDD1A969B1CB9A249161958903A4ADB5'
    }
    # iterate through each PMID and generate the url with the specific PMID
    for pmid in pmids:
        print("Extracting Annotations for PMID: {0}".format(pmid))
        annotations_url = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=MED%3A" + pmid + "&type=Gene%20Function&section=Abstract&provider=HES-SO_SIB&format=JSON"
        # send HTTP get request and append the JSON formatted output to the list
        response = requests.request("GET", annotations_url, headers=headers, data=payload, verify=False)
        annotations.append(response.json())
        # give a delay of 0.5s to avoid server from crashing due to continuous request to the API server
        print("Done")
        time.sleep(0.5)
    # finally, return the responses
    return annotations


def get_go_terms(pmids):
    """
    This method extracts the GO terms for each PMID using the API
    :param pmids: the list of PMIDs
    :return: the list containing the HTTP response with GO terms for each PMID
    """
    # initialize an empty list to hold the results
    go_terms = []
    # for this API, payload is empty
    payload = {}
    headers = {
        'Cookie': 'JSESSIONID=EDD1A969B1CB9A249161958903A4ADB5'
    }
    # iterate through each PMID
    for pmid in pmids:
        print("Extracting Go Terms for PMID: {0}".format(pmid))
        # generate the url to get the GO terms for the specific PMID
        annotations_url = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=MED%3A" + pmid + "&type=Gene%20Ontology&section=Abstract&provider=Europe%20PMC&format=JSON"
        # send HTTP Get request
        response = requests.request("GET", annotations_url, headers=headers, data=payload, timeout=2, verify=False)
        # append to the results list, the JSON formatted output
        go_terms.append(response.json())
        print("Done.")
        # give a delay of 0.5s to avoid server from getting crashed due to continuous requests to API server
        time.sleep(0.5)
    # finally, return the response
    return go_terms


def get_consolidated_annotations(pmids, annotations, go_terms, json_data):
    """
    This function consolidates the annotations based on whether there's annotation for a
    PMID or not. If result is present from the annotations API for a particular PMID, that
    result is used. Otherwise, the corresponding Go term is used for that PMID
    :param pmids: the list of PMIDs
    :param annotations: the list of annotation responses
    :param go_terms: the list of Go terms responses
    :return: a consolidated dictionary with PMID as key and annotations as its value
    """
    # initialize an empty dictionary to hold the results
    pmid_annotations = {}
    # iterate through each of the annotations
    for i in range(0, len(pmids)):
        # extract the PMID, this is available in the 'extId' attribute
        if len(annotations[i]) != 0:
            pmid = annotations[i][0]['extId']
            # extract the annotation
            annotation = annotations[i][0]['annotations']
            # if annotation is empty for that specific PMID,
            if len(annotation) == 0:
                # get the go term annotation for the specific PMID and
                # add this to the key-value pair in the dictionary
                go_trm = go_terms[i][0]['annotations']
                updated_go_terms = []
                for term in go_trm:
                    if term['exact'] != '':
                        exact_term = term['exact']
                        # extract the sentence index
                        sentence_data = json.loads(json_data[i])['text']
                        sentences = sentence_data.split('.')
                        response = [(i, sentence) for i, sentence in enumerate(sentences)
                                     if exact_term in sentence]
                        if len(response) != 0:
                            start_index, sentence = response[0]
                            term['sentenceIndex'] = start_index
                            # extract charStart
                            term['charStart'] = sentence_data.find(exact_term)
                            # extract charEnd
                            term['charEnd'] = term['charStart'] + len(exact_term)
                        updated_go_terms.append(term)

                pmid_annotations[pmid] = updated_go_terms
            else:
                # if annotation is present, add this to the key-value pair
                annotation_list = []
                for annot in annotation:
                    # this will be a dict
                    if annot['exact'] != "":
                        # get the required params
                        # sentence index
                        exact_term = annot['exact']
                        sentence_data = json.loads(json_data[i])['text']
                        sentences = sentence_data.split(".")
                        exact_term = exact_term.replace(".", "")
                        start_index, sentence = [(idx, s) for idx, s in enumerate(sentences) if exact_term in s][0]
                        annot['sentenceIndex'] = start_index
                        annot['charStart'] = sentence_data.find(exact_term)
                        annot['charEnd'] = annot['charStart'] + len(exact_term)
                        annotation_list.append(annot)
                pmid_annotations[pmid] = annotation_list
    # finally, return the result
    return pmid_annotations


def generate_result(pmids, pgenn_results, pmid_annotations, json_data):
    global date
    """
    This function will generate a .csv file with all the extracted results from the APIs
    :return: None
    """
    # these are the columns for the resulting CSV file
    columns = ['PMID', 'PGennGeneName', 'UniprotID', 'EntrezID', 'Annotations',
               'Title', 'Text', 'Sentence', 'Mesh', 'Type']
    # initialize empty list to hold the consolidated entries
    records = []

    # iterate through each PMID
    for idx, pmid in enumerate(pmids):
        # create an empty dictionary holding entry information
        entry = {}
        # add to the dictionary, the PMID
        entry['PMID'] = pmid
        # iterate through each pgenn_result key-value pair
        if pmid in pgenn_results.keys():
            # set the PGennGeneName and UniprotID to be empty initially
            entry['PGennGeneName'] = []
            # iterate through each results for the particular PMID
            genes_info = []
            for item in pgenn_results[pmid]:
                # add the result as comma separated values for PGennGeneName and UniprotID
                gene_name_details = {}
                gene_name = item[0]
                gene_name_details['Name'] = gene_name
                gene_name_details['UniprotID'] = item[1]
                gene_name_details['EntrezID'] = item[2]
                # extract the sentence index
                sentence_data = json.loads(json_data[idx])['text']
                sentences = sentence_data.split('.')
                response = [(i, sentence) for i, sentence in enumerate(sentences)
                                         if gene_name in sentence]
                if len(response) != 0:
                    start_index, sentence = response[0]
                    gene_name_details['sentenceIndex'] = start_index
                    # extract charStart
                    gene_name_details['charStart'] = sentence_data.find(gene_name)
                    # extract charEnd
                    gene_name_details['charEnd'] = gene_name_details['charStart'] + len(gene_name)
                genes_info.append(gene_name_details)
            entry['PGennGeneName'] = genes_info
        else:
            entry['PGennGeneName'] = []

        # add annotation info
        if pmid in pmid_annotations.keys():
            entry['Annotations'] = pmid_annotations[pmid]
        else:
            entry['Annotations'] = ""

        # append the json input data also in the records
        json_entry = json.loads(json_data[idx])
        try:
            entry['Title'] = json_entry['title']
            entry['Text'] = json_entry['text']
            entry['Sentence'] = json_entry['sentence']
            entry['Mesh'] = json_entry['mesh']
            entry['Type'] = json_entry['type']
        except Exception as ex:
            print(f"Error: {ex}")

        # append the entry to the list
        if len(entry['PGennGeneName']) != 0 and len(entry['Annotations']) != 0:
            records.append(entry)

    # create a pandas DataFrame with the list of dictionaries with entries
    df = pd.DataFrame(data=records, columns=columns)

    # save the DataFrame as csv File
    df.to_csv(f'results_{date}.csv', mode='a', header=False)
    # return the list of dictionary to save in database
    generate_database.save(records)


def extract(json_data):
    """
    This is the main function that invokes the helper functions
    :return: None
    """
    try:
        # invoke the function to get the list of PMIDs
        print("Extracting PMIDs...")
        # the below line extracts PMIDs from JSON file
        pmids = extract_pmids_from_json(json_data)
        # uncomment below line if API should be Pubmed
        # pmids = extract_pmids()
        print("Completed.")
        # use this PMID list to extract PGenn API responses
        print("Extracting PGenn Responses...")
        pgenn_responses = extract_pgenn_response(pmids)
        print("Completed.")
        # parse this response to generate a dictionary of gene name and UniprotIDs and Entrez IDs for the PMIDs
        print("Parsing Gene Name, UniprotID and EntrezID...")
        gene_uniprot_entrez_map = get_gene_uniprot(pgenn_responses)
        print("Completed.")
        # # use the PMID list to extract Pubtator API responses
        # print("Extracting Pubtator Responses...")
        # pubtator_responses = extract_pubtator_response(pmids)
        # print("Completed.")
        # # parse this response to generate a dictionary of gene name and EntrezIDs for the PMIDs
        # print("Parsing Gene Name and EntrezID...")
        # gene_entrez_id_map = get_gene_entrez_id(pubtator_responses)
        # print("Completed.")
        # get the annotations for all PMIDs
        print("Extracting Annotations...")
        annotations = get_annotations(pmids)
        print("Completed.")
        # get the Go terms for all PMIDs
        print("Extracting Go Terms...")
        go_terms = get_go_terms(pmids)
        print("Completed.")
        # use the annotations and go terms to consolidate the annotations for each PMID
        print("Generating Consolidated Annotations...")
        consolidated_annotations = get_consolidated_annotations(pmids, annotations, go_terms, json_data)
        print("Completed.")
        # use all the results to generate the consolidated file
        print("Generating result...")
        generate_result(pmids, gene_uniprot_entrez_map, consolidated_annotations, json_data)
        print("Completed.")
    except Exception as ex:
        print("Error: {0}".format(ex))


if __name__ == '__main__':
    parts = ['{"docId": "30309330", "title": {"charStart": 0, "charEnd": 87}, "text": "CsINV5, a tea vacuolar invertase gene enhances cold tolerance in transgenic Arabidopsis. BACKGROUND: Vacuolar invertases (VINs) have been reported to regulate plant growth and development and respond to abiotic stresses such as drought and cold. With our best knowledge, the functions of VIN genes little have been reported in tea plant (Camellia sinensis L.). Therefore, it is necessary to develop research in this field. RESULTS: Here, we identified a VIN gene, CsINV5, which was induced by cold acclimation and sugar treatments in the tea plant. Histochemical assays results showed that the 1154\\u00a0bp 5\'-flanking sequence of CsINV5 drove \\u03b2-glucuronidase (GUS) gene expression in roots, stems, leaves, flowers and siliques of transgenic Arabidopsis during different developmental stages. Moreover, promoter deletion analysis results revealed that an LTRE-related motif (CCGAAA) and a WBOXHVISO1 motif (TGACT) within the promoter region of CsINV5 were the core cis-elements in response to low temperature and sugar signaling, respectively. In addition, overexpression of CsINV5 in Arabidopsis promoted taproot and lateral root elongation through glucose-mediated effects on auxin signaling. Based on physiological and RNA-seq analysis, we found that overexpression of CsINV5 improved cold tolerance in transgenic Arabidopsis mainly by increasing the contents of glucose and fructose, the corresponding ratio of hexose to sucrose, and the transcription of osmotic-stress-related genes (P5CS1, P5CS2, AtLEA3, COR413-PM1 and COR15B) to adjust its osmotic potential. CONCLUSIONS: Comprehensive experimental results suggest that overexpression of CsINV5 may enhance the cold tolerance of plant through the modification of cellular sugar compounds contents and osmotic regulation related pathways.", "sentence": [{"charStart": 0, "charEnd": 87, "index": 0}, {"index": 1, "charEnd": 244, "charStart": 89}, {"index": 2, "charEnd": 359, "charStart": 246}, {"index": 3, "charEnd": 421, "charStart": 361}, {"index": 4, "charEnd": 547, "charStart": 423}, {"index": 5, "charEnd": 786, "charStart": 549}, {"index": 6, "charEnd": 1037, "charStart": 788}, {"index": 7, "charEnd": 1188, "charStart": 1039}, {"index": 8, "charEnd": 1560, "charStart": 1190}, {"index": 9, "charEnd": 1789, "charStart": 1562}], "mesh": ["Arabidopsis", "Cold Temperature", "Gene Expression Regulation, Plant", "Plant Proteins", "Plants, Genetically Modified", "Tea", "beta-Fructofuranosidase"], "type": ["Journal Article"]}', '{"docId": "30312933", "title": {"charStart": 0, "charEnd": 169}, "text": "In silico characterization and transcriptional modulation of phenylalanine ammonia lyase (PAL) by abiotic stresses in the medicinal orchid Vanda coerulea Griff. ex Lindl. Phenylalanine ammonia lyase (PAL) is the first enzyme of phenylpropanoid pathway. In the present study, a full-length PAL transcript from Vanda coerulea Griff. ex Lindl. (Family: Orchidaceae) was isolated and characterized. It was found that complete PAL transcript of V. coerulea (VcPAL; Gene Bank no. MG745168) contained 2175 bp with the open reading frame (ORF) of 2112 bp, encoding 703 amino acid residues. The multiple sequence alignment showed that VcPAL protein had 81% identity with that of the orchid, Bromheadia finlaysoniana. Phylogenetic analysis also disclosed that VcPAL shared the same evolutionary relationship with PAL proteins of other orchid species and to be closely related to that of other angiosperm species as well. The three-dimensional structure of VcPAL was found to be homo-tetrameric in nature consisting of four identical subunits with a molecular mass of 75\\u202fkDa per subunit. In silico characterization revealed the deduced protein to be a stable protein, comprising three major functional domains as reported in PAL proteins of other species. The transcription profiling of VcPAL exhibited the highest expression level to be present in the in vitro - raised leaf and root samples as compared to that of the ex vitro plant. The differential expression of VcPAL transcript was observed to be up-regulated by different types of abiotic stresses like wounding, cold, UV-B, salinity, and down-regulated by dark treatment. The study also exhibited that the VcPAL enzyme activity was directly proportional to the gene expression after the tissues were subjected to salinity and wounding stresses wherein a 1.7- fold increase in the enzyme activity was recorded in the leaf tissues exposed to salinity stress. A positive correlation could be found between the enzyme activity and the accumulation of phenylpropanoids such as total phenolic and flavonoid contents with R2\\u202f=\\u202f0.85 and 0.842 respectively.", "sentence": [{"charStart": 0, "charEnd": 159, "index": 0}, {"index": 1, "charEnd": 169, "charStart": 161}, {"index": 2, "charEnd": 251, "charStart": 171}, {"index": 3, "charEnd": 329, "charStart": 253}, {"index": 4, "charEnd": 339, "charStart": 331}, {"index": 5, "charEnd": 393, "charStart": 341}, {"index": 6, "charEnd": 472, "charStart": 395}, {"index": 7, "charEnd": 580, "charStart": 474}, {"index": 8, "charEnd": 706, "charStart": 582}, {"index": 9, "charEnd": 909, "charStart": 708}, {"index": 10, "charEnd": 1075, "charStart": 911}, {"index": 11, "charEnd": 1243, "charStart": 1077}, {"index": 12, "charEnd": 1423, "charStart": 1245}, {"index": 13, "charEnd": 1617, "charStart": 1425}, {"index": 14, "charEnd": 1902, "charStart": 1619}, {"index": 15, "charEnd": 2094, "charStart": 1904}], "mesh": ["Amino Acid Sequence", "Computer Simulation", "Molecular Conformation", "Orchidaceae", "Phenylalanine Ammonia-Lyase", "Sequence Alignment", "Transcription, Genetic"], "type": ["Journal Article"]}', '{"docId": "30312954", "title": {"charStart": 0, "charEnd": 95}, "text": "Overexpression of Brassica campestris BcICE1 gene increases abiotic stress tolerance in tobacco. In this study, a cDNA of ICE1 (inducer of CBF expression 1) gene, named BcICE1, was isolated from Brassica campestris \'Longyou 6\'. The deduced protein has 499 amino acids containing a typical bHLH domain and is highly identical with AtICE1 (85.9%) from Arabidopsis thaliana. BcICE1 is located in the nucleus. The activities of SOD, CAT, POD, and APX and the transcriptional levels of SOD, CAT, and POD genes were higher in BcICE1-transgenic tobacco than in wild-type (WT) tobacco under cold stress. Compared with WT tobacco, proline, soluble sugar, and chlorophyll were enhanced, whereas malondialdehyde and relative conductivity were decreased in BcICE1-transgenic tobacco. The overexpression of BcICE1 in tobacco increased the expression of CBF1, CBF2, and other stress-related genes. Moreover, under salt and PEG (25%) stress, the activities of APX and GPX and content of soluble sugar and chlorophyll in BcICE1-transgenic tobacco were higher than those in WT tobacco. Our results suggest that BcICE1 plays an important role in the response to abiotic stress.", "sentence": [{"charStart": 0, "charEnd": 95, "index": 0}, {"index": 1, "charEnd": 226, "charStart": 97}, {"index": 2, "charEnd": 370, "charStart": 228}, {"index": 3, "charEnd": 404, "charStart": 372}, {"index": 4, "charEnd": 594, "charStart": 406}, {"index": 5, "charEnd": 770, "charStart": 596}, {"index": 6, "charEnd": 882, "charStart": 772}, {"index": 7, "charEnd": 1067, "charStart": 884}, {"index": 8, "charEnd": 1158, "charStart": 1069}], "mesh": ["Adaptation, Physiological", "Amino Acid Sequence", "Antioxidants", "Brassica", "Cell Nucleus", "Chlorophyll", "Cloning, Molecular", "Gene Expression Regulation, Plant", "Genes, Plant", "Hydrogen Peroxide", "Malondialdehyde", "Plant Proteins", "Plants, Genetically Modified", "Proline", "Stress, Physiological", "Subcellular Fractions", "Tobacco"], "type": ["Journal Article"]}', '{"docId": "30312968", "title": {"charStart": 0, "charEnd": 187}, "text": "A stress-associated protein, LmSAP, from the halophyte Lobularia maritima provides tolerance to heavy metals in tobacco through increased ROS scavenging and metal detoxification processes. Agricultural soil pollution by heavy metals is a severe global ecological problem. We recently showed that overexpression of LmSAP, a member of the stress-associated protein (SAP) gene family isolated from Lobularia maritima, in transgenic tobacco led to enhanced tolerance to abiotic stress. In this study, we characterised the response of LmSAP transgenic tobacco plants to metal stresses (cadmium (Cd), copper (Cu), manganese (Mn), and zinc (Zn)). In L. maritima, LmSAP expression increased after 12\\u2009h of treatment with these metals, suggesting its involvement in the plant response to heavy metal stress. LmSAP transgenic tobacco plants subjected to these stress conditions were healthy, experienced higher seedling survival rates, and had longer roots than non-transgenic plants (NT). However, they exhibited higher tolerance towards cadmium and manganese than towards copper and zinc. LmSAP-overexpressing tobacco seedlings accumulated more cadmium, copper, and manganese compared with NT plants, but displayed markedly decreased hydrogen peroxide (H2O2) and lipid peroxidation levels after metal treatment. Activities of the antioxidant enzymes superoxide dismutase (SOD), catalase (CAT), and peroxidase (POD) were significantly higher in transgenic plants than in NT plants after exposure to metal stress. LmSAP overexpression also enhanced the transcription of several genes encoding metallothioneins (Met1, Met2, Met3, Met4, and Met5), a copper transport protein CCH, a Cys and His-rich domain-containing protein RAR1 (Rar1), and a ubiquitin-like protein 5 (PUB1), which are involved in metal tolerance in tobacco. Our findings indicate that LmSAP overexpression in tobacco enhanced tolerance to heavy metal stress by protecting the plant cells against oxidative stress, scavenging reactive oxygen species (ROS), and decreasing the intracellular concentration of free heavy metals through its effect on metal-binding proteins in the cytosol.", "sentence": [{"charStart": 0, "charEnd": 187, "index": 0}, {"index": 1, "charEnd": 270, "charStart": 189}, {"index": 2, "charEnd": 480, "charStart": 272}, {"index": 3, "charEnd": 638, "charStart": 482}, {"index": 4, "charEnd": 796, "charStart": 640}, {"index": 5, "charEnd": 977, "charStart": 798}, {"index": 6, "charEnd": 1078, "charStart": 979}, {"index": 7, "charEnd": 1301, "charStart": 1080}, {"index": 8, "charEnd": 1501, "charStart": 1303}, {"index": 9, "charEnd": 1812, "charStart": 1503}, {"index": 10, "charEnd": 2139, "charStart": 1814}], "mesh": ["Brassicaceae", "Genes, Plant", "Metals, Heavy", "Oxidative Stress", "Plant Proteins", "Plants, Genetically Modified", "Reactive Oxygen Species", "Real-Time Polymerase Chain Reaction", "Salt-Tolerant Plants", "Soil Pollutants", "Tobacco"], "type": ["Journal Article"]}', '{"docId": "30314273", "title": {"charStart": 0, "charEnd": 100}, "text": "Overexpression of a Novel ROP Gene from the Banana (MaROP5g) Confers Increased Salt Stress Tolerance. Rho-like GTPases from plants (ROPs) are plant-specific molecular switches that are crucial for plant survival when subjected to abiotic stress. We identified and characterized 17 novel ROP proteins from Musa acuminata (MaROPs) using genomic techniques. The identified MaROPs fell into three of the four previously described ROP groups (Groups II\\u207bIV), with MaROPs in each group having similar genetic structures and conserved motifs. Our transcriptomic analysis showed that the two banana genotypes tested, Fen Jiao and BaXi Jiao, had similar responses to abiotic stress: Six genes (MaROP-3b, -5a, -5c, -5f, -5g, and -6) were highly expressed in response to cold, salt, and drought stress conditions in both genotypes. Of these, MaROP5g was most highly expressed in response to salt stress. Co-localization experiments showed that the MaROP5g protein was localized at the plasma membrane. When subjected to salt stress, transgenic Arabidopsis thaliana overexpressing MaROP5g had longer primary roots and increased survival rates compared to wild-type A. thaliana. The increased salt tolerance conferred by MaROP5g might be related to reduced membrane injury and the increased cytosolic K\\u207a/Na\\u207a ratio and Ca2+ concentration in the transgenic plants as compared to wild-type. The increased expression of salt overly sensitive (SOS)-pathway genes and calcium-signaling pathway genes in MaROP5g-overexpressing A. thaliana reflected the enhanced tolerance to salt stress by the transgenic lines in comparison to wild-type. Collectively, our results suggested that abiotic stress tolerance in banana plants might be regulated by multiple MaROPs, and that MaROP5g might enhance salt tolerance by increasing root length, improving membrane injury and ion distribution.", "sentence": [{"charStart": 0, "charEnd": 100, "index": 0}, {"index": 1, "charEnd": 244, "charStart": 102}, {"index": 2, "charEnd": 353, "charStart": 246}, {"index": 3, "charEnd": 533, "charStart": 355}, {"index": 4, "charEnd": 818, "charStart": 535}, {"index": 5, "charEnd": 890, "charStart": 820}, {"index": 6, "charEnd": 988, "charStart": 892}, {"index": 7, "charEnd": 1163, "charStart": 990}, {"index": 8, "charEnd": 1372, "charStart": 1165}, {"index": 9, "charEnd": 1616, "charStart": 1374}, {"index": 10, "charEnd": 1859, "charStart": 1618}], "mesh": ["Adaptation, Biological", "Biomarkers", "Computational Biology", "Conserved Sequence", "Gene Expression Regulation, Plant", "Multigene Family", "Musa", "Nucleotide Motifs", "Phenotype", "Phylogeny", "Plants, Genetically Modified", "Reactive Oxygen Species", "Reproducibility of Results", "Salt Stress", "Salt Tolerance", "Signal Transduction", "Stress, Physiological", "rho GTP-Binding Proteins"], "type": ["Journal Article"]}', '{"docId": "30315246", "title": {"charStart": 0, "charEnd": 130}, "text": "Comparative transcriptomic analysis reveals common molecular factors responsive to heat and drought stress in Agrostis stolonifera. Heat and drought stress are primary abiotic stresses confining growth of cool-season grass species during summer. The objective of this study was to identify common molecular factors and metabolic pathways associated with heat and drought responses in creeping bentgrass (Agrostis stolonifera) by comparative analysis of transcriptomic profiles between plants exposed to heat and drought stress. Plants were exposed to heat stress (35/30\\u2009\\u00b0C day/night temperature) or drought stress by withholding irrigation for 21 d in growth chambers. Transcriptomic profiling by RNA-seq in A. stolonifera (cv. \'Penncross\') found 670 commonly up-regulated and 812 commonly down-regulated genes by heat and drought stress. Transcriptional up-regulations of differentially expressed genes (DEGs) due to heat and drought stress include genes that were highly enriched in oxylipin biosynthetic process and proline biosynthetic process. Transcriptional down-regulations of genes under heat and drought stress were highly enriched and involved in thiamine metabolic process and calcium sensing receptor. These commonly-regulated genes by heat and drought stress identified in A. stolonifera suggested that drought and heat responses shared such common molecular factors and pathways, which could be potential candidate genes for genetic modification of improving plant tolerance to the combined heat and drought stress.", "sentence": [{"charStart": 0, "charEnd": 130, "index": 0}, {"index": 1, "charEnd": 244, "charStart": 132}, {"index": 2, "charEnd": 526, "charStart": 246}, {"index": 3, "charEnd": 667, "charStart": 528}, {"index": 4, "charEnd": 726, "charStart": 669}, {"index": 5, "charEnd": 837, "charStart": 728}, {"index": 6, "charEnd": 1047, "charStart": 839}, {"index": 7, "charEnd": 1213, "charStart": 1049}, {"index": 8, "charEnd": 1529, "charStart": 1215}], "mesh": ["Agrostis", "Down-Regulation", "Droughts", "Electrolytes", "Gene Expression Profiling", "Gene Expression Regulation, Plant", "Gene Ontology", "Hot Temperature", "Plant Leaves", "Reproducibility of Results", "Stress, Physiological", "Up-Regulation", "Water"], "type": ["Comparative Study", "Journal Article", "Research Support, Non-U.S. Gov\'t"]}', '{"docId": "30316294", "title": {"charStart": 0, "charEnd": 118}, "text": "A Glycine soja group S2 bZIP transcription factor GsbZIP67 conferred bicarbonate alkaline tolerance in Medicago sativa. BACKGROUND: Even though bicarbonate alkaline stress is a serious threat to crop growth and yields, it attracts much fewer researches than high salinity stress. The basic leucine zipper (bZIP) transcription factors have been well demonstrated to function in diverse abiotic stresses; however, their biological role in alkaline tolerance still remains elusive. In this study, we functionally characterized a bZIP gene from Glycine soja GsbZIP67 in bicarbonate alkaline stress responses. RESULTS: GsbZIP67 was initially identified as a putative bicarbonate responsive gene, on the basis of previous RNA-seq data of 50\\u00a0mM NaHCO3-treated Glycine soja roots. GsbZIP67 protein possessed a conserved bZIP domain, and belonged to the group S2 bZIP, which is yet less well-studied. Our studies showed that GsbZIP67 targeted to nucleus in Arabidopsis protoplasts, and displayed transcriptional activation activity in yeast cells. The quantitative real-time PCR analyses unraveled the bicarbonate stress responsive expression and tissue specific expression of GsbZIP67 in wild soybean. Further phenotypic analysis illustrated that GsbZIP67 overexpression in alfalfa promoted plant growth under bicarbonate alkaline stress, as evidenced by longer roots and shoots. Furthermore, GsbZIP67 overexpression also modified the physiological indices of transgenic alfalfa under bicarbonate alkaline stress. In addition, the expression levels of several stress responsive genes were also augmented by GsbZIP67 overexpression. CONCLUSIONS: Collectively, in this study, we demonstrated that GsbZIP67 acted as a positive regulator of plant tolerance to bicarbonate alkaline stress. These results provide direct genetic evidence of group S2 bZIPs in bicarbonate alkaline stress, and will facilitate further studies concerning the cis-elements and/or downstream genes targeted by GsbZIP67 in stress responses.", "sentence": [{"charStart": 0, "charEnd": 118, "index": 0}, {"index": 1, "charEnd": 278, "charStart": 120}, {"index": 2, "charEnd": 477, "charStart": 280}, {"index": 3, "charEnd": 603, "charStart": 479}, {"index": 4, "charEnd": 771, "charStart": 605}, {"index": 5, "charEnd": 890, "charStart": 773}, {"index": 6, "charEnd": 1037, "charStart": 892}, {"index": 7, "charEnd": 1192, "charStart": 1039}, {"index": 8, "charEnd": 1370, "charStart": 1194}, {"index": 9, "charEnd": 1504, "charStart": 1372}, {"index": 10, "charEnd": 1622, "charStart": 1506}, {"index": 11, "charEnd": 1775, "charStart": 1624}, {"index": 12, "charEnd": 2001, "charStart": 1777}], "mesh": ["Alkalies", "Amino Acid Sequence", "Basic-Leucine Zipper Transcription Factors", "Bicarbonates", "Cell Nucleus", "Gene Expression Regulation, Plant", "Genes, Reporter", "Medicago sativa", "Phenotype", "Phylogeny", "Plant Proteins", "Plant Roots", "Plant Shoots", "Plants, Genetically Modified", "Protein Transport", "Sequence Alignment", "Soybeans", "Stress, Physiological"], "type": ["Journal Article"]}', '{"docId": "30324736", "title": {"charStart": 0, "charEnd": 64}, "text": "ABA signalling is regulated by the circadian clock component LHY.", "sentence": [{"charStart": 0, "charEnd": 64, "index": 0}], "mesh": ["Abscisic Acid", "Arabidopsis", "Arabidopsis Proteins", "Circadian Clocks", "Signal Transduction"], "type": ["Journal Article", "Comment"]}', '{"docId": "30326419", "title": {"charStart": 0, "charEnd": 88}, "text": "Effect of heat root stress and high salinity on glucosinolates metabolism in wild rocket. Wild rocket (Diplotaxis tenuifolia L.) is a leafy vegetable appreciated for its characteristic sensory properties which are mainly due to the presence of glucosinolates (GSLs). Short-term exposure to abiotic stresses can induce physiological responses and transcriptional changes which involve GSLs. For this reason, the aim of this work was to study the mechanisms of regulation of GSLs metabolism in rocket subjected to heat stress (40\\u2009\\u00b0C) and high salinity (200\\u2009mM NaCl) imposed for up to 48\\u2009h. GSLs levels and the expression of methylthioalkylmalate synthase1 (DtMAM1), cytochromeP79F1 (DtCYP79F1), cytochromeP45083A1 (DtCYP83A1), cytosolic-sulfotransferase5b (DtST5b), cytosolic-sulfotransferase5c (DtST5c), flavinmono-oxygenase (DtFMO), myrosinase (DtMyro) and thio-methyl transferase (DtTMT) were analyzed under stress conditions. In addition, the effect on chlorophyll and glucose levels, as well as on chlorophyll a fluorescence were evaluated. Chlorophyll and chlorophyll fluorescence were not affected by the short-term application of stresses. Glucose levels in roots were doubled in response to high salinity, while, in the same organ, GSLs were three fold lower in response to both stresses. The relative content of several aliphatic GSLs was significantly reduced in leaves as a response to both stresses. A key role in GSLs metabolism and in the response to salinity is hypothesized for the gene DtTMT, as it showed an increment in transcripts accumulation (three-fold) consistent with the decrement in the GSLs levels found in salt-exposed leaves and roots. The results obtained in this study can be used in breeding programmes aiming to enhance rocket sensory quality and to improve the resistance to abiotic stresses.", "sentence": [{"charStart": 0, "charEnd": 88, "index": 0}, {"index": 1, "charEnd": 265, "charStart": 90}, {"index": 2, "charEnd": 388, "charStart": 267}, {"index": 3, "charEnd": 926, "charStart": 390}, {"index": 4, "charEnd": 1042, "charStart": 928}, {"index": 5, "charEnd": 1144, "charStart": 1044}, {"index": 6, "charEnd": 1294, "charStart": 1146}, {"index": 7, "charEnd": 1409, "charStart": 1296}, {"index": 8, "charEnd": 1663, "charStart": 1411}, {"index": 9, "charEnd": 1825, "charStart": 1665}], "mesh": ["Brassicaceae", "Chlorophyll", "Chromatography, High Pressure Liquid", "Gene Expression Regulation, Plant", "Glucose", "Glucosinolates", "Hot Temperature", "Plant Leaves", "Plant Roots", "Salt Stress", "Stress, Physiological"], "type": ["Journal Article"]}', '{"docId": "30332983", "title": {"charStart": 0, "charEnd": 94}, "text": "Genome-wide identification and analysis of the COI gene family in wheat (Triticum aestivum L.). BACKGROUND: COI (CORONATINE INSENSITIVE), an F-box component of the Skp1-Cullin-F-box protein (SCFCOI1) ubiquitin E3 ligase, plays important roles in the regulation of plant growth and development. Recent studies have shown that COIs are involved in pollen fertility. In this study, we identified and characterized COI genes in the wheat genome and analyzed expression patterns under abiotic stress. RESULTS: A total of 18 COI candidate sequences for 8 members of COI gene family were isolated in wheat (Triticum aestivum L.). Phylogenetic and structural analyses showed that these COI genes could be divided into seven distinct subfamilies. The COI genes showed high expression in stamens and glumes. The qRT-PCR results revealed that wheat COIs were involved in several abiotic stress responses and anther/glume dehiscence in the photoperiod-temperature sensitive genic male sterile (PTGMS) wheat line BS366. CONCLUSIONS: The structural characteristics and expression patterns of the COI gene family in wheat as well as the stress-responsive and differential tissue-specific expression profiles of each TaCOI gene were examined in PTGMS wheat line BS366. In addition, we examined SA- and MeJA-induced gene expression in the wheat anther and glume to investigate the role of COI in the JA signaling pathway, involved in the regulation of abnormal anther dehiscence in the PTGMS wheat line. The results of this study contribute novel and detailed information about the TaCOI gene family in wheat and could be used as a benchmark for future studies of the molecular mechanisms of PTGMS in other crops.", "sentence": [{"charStart": 0, "charEnd": 94, "index": 0}, {"index": 1, "charEnd": 292, "charStart": 96}, {"index": 2, "charEnd": 362, "charStart": 294}, {"index": 3, "charEnd": 494, "charStart": 364}, {"index": 4, "charEnd": 621, "charStart": 496}, {"index": 5, "charEnd": 736, "charStart": 623}, {"index": 6, "charEnd": 796, "charStart": 738}, {"index": 7, "charEnd": 1005, "charStart": 798}, {"index": 8, "charEnd": 1251, "charStart": 1007}, {"index": 9, "charEnd": 1485, "charStart": 1253}, {"index": 10, "charEnd": 1695, "charStart": 1487}], "mesh": ["Cyclopentanes", "Gene Expression Profiling", "Genome, Plant", "Genomics", "Organ Specificity", "Oxylipins", "Phylogeny", "Promoter Regions, Genetic", "Signal Transduction", "Triticum", "Ubiquitin-Protein Ligases"], "type": ["Journal Article"]}', '{"docId": "30333149", "title": {"charStart": 0, "charEnd": 111}, "text": "The SUMO E3 Ligase MdSIZ1 Targets MdbHLH104 to Regulate Plasma Membrane H+-ATPase Activity and Iron Homeostasis. SIZ1 (a SIZ/PIAS-type SUMO E3 ligase)-mediated small ubiquitin-like modifier (SUMO) modification of target proteins is important for various biological processes related to abiotic stress resistance in plants; however, little is known about its role in resistance toward iron (Fe) deficiency. Here, the SUMO E3 ligase MdSIZ1 was shown to be involved in the plasma membrane (PM) H+-ATPase-mediated response to Fe deficiency. Subsequently, a basic helix-loop-helix transcription factor, MdbHLH104 (a homolog of Arabidopsis bHLH104 in apple), which acts as a key component in regulating PM H+-ATPase-mediated rhizosphere acidification and Fe uptake in apples (Malus domestica), was identified as a direct target of MdSIZ1. MdSIZ1 directly sumoylated MdbHLH104 both in vitro and in vivo, especially under conditions of Fe deficiency, and this sumoylation was required for MdbHLH104 protein stability. Double substitution of K139R and K153R in MdbHLH104 blocked MdSIZ1-mediated sumoylation in vitro and in vivo, indicating that the K139 and K153 residues were the principal sites of SUMO conjugation. Moreover, the transcript level of the MdSIZ1 gene was substantially induced following Fe deficiency. MdSIZ1 overexpression exerted a positive influence on PM H+-ATPase-mediated rhizosphere acidification and Fe uptake. Our findings reveal an important role for sumoylation in the regulation of PM H+-ATPase-mediated rhizosphere acidification and Fe uptake during Fe deficiency in plants.", "sentence": [{"charStart": 0, "charEnd": 111, "index": 0}, {"index": 1, "charEnd": 404, "charStart": 113}, {"index": 2, "charEnd": 535, "charStart": 406}, {"index": 3, "charEnd": 831, "charStart": 537}, {"index": 4, "charEnd": 1008, "charStart": 833}, {"index": 5, "charEnd": 1207, "charStart": 1010}, {"index": 6, "charEnd": 1308, "charStart": 1209}, {"index": 7, "charEnd": 1425, "charStart": 1310}, {"index": 8, "charEnd": 1594, "charStart": 1427}], "mesh": ["Cell Membrane", "Iron", "Malus", "Proton-Translocating ATPases", "RNA, Messenger", "Rhizosphere", "Sumoylation", "Ubiquitins"], "type": ["Journal Article", "Research Support, Non-U.S. Gov\'t"]}', '{"docId": "30335565", "title": {"charStart": 0, "charEnd": 81}, "text": "Protein kinase CK2\\u03b1 subunits constitutively activate ABA signaling in Arabidopsis. Protein kinase CK2 (formerly known as casein kinase II), a Ser/Thr protein kinase highly conserved in eukaryotes, is essential for cell survival by regulating\\u00a0a wide range of plant growth, development, and stress responses. A growing body of evidence has shown a link between CK2 and abscisic acid (ABA) signaling in response to abiotic stress. However, the roles of CK2 subunits in ABA signaling remain unclear in plants. Our recent work in Arabidopsis thaliana has revealed that CK2\\u03b1 and CK2\\u03b2 subunits inversely modulate ABA signal output. Here, we examine the roles of CK2\\u03b1s, by assessing how CK2\\u03b1s affect ABA signaling. Together with the previous findings, our mutant and transient expression analyses demonstrate that CK2\\u03b1s positively modulate ABA signaling through the core ABA signaling pathway in the presence of ABA, though the positive effect of CK2\\u03b1s are much smaller than that of core ABA signaling components in ABA response. In addtion, our current and previous findings also suggest that CK2\\u03b1s play a role in maintaining constitutively active ABA signaling even in the absence of ABA independently of the core ABA signaling pathway. Thus, we found that CK2\\u03b1s constitutively activate ABA signaling in the presence or absence of ABA in a different manner in Arabidopsis plants.", "sentence": [{"charStart": 0, "charEnd": 81, "index": 0}, {"index": 1, "charEnd": 305, "charStart": 83}, {"index": 2, "charEnd": 426, "charStart": 307}, {"index": 3, "charEnd": 504, "charStart": 428}, {"index": 4, "charEnd": 623, "charStart": 506}, {"index": 5, "charEnd": 705, "charStart": 625}, {"index": 6, "charEnd": 1020, "charStart": 707}, {"index": 7, "charEnd": 1229, "charStart": 1022}, {"index": 8, "charEnd": 1372, "charStart": 1231}], "mesh": ["Abscisic Acid", "Arabidopsis", "Arabidopsis Proteins", "Casein Kinase II", "Gene Expression Regulation, Plant", "Signal Transduction", "Stress, Physiological"], "type": ["Journal Article", "Research Support, Non-U.S. Gov\'t"]}', '{"docId": "30340051", "title": {"charStart": 0, "charEnd": 132}, "text": "Cloning and functional characterization of the Na+/H+ antiporter (NHX1) gene promoter from an extreme halophyte Salicornia brachiata. Salinity is one of the major abiotic stresses which affect plant growth and productivity by imposing dual stress, ionic and osmotic stress, on plants. Halophytes which are adapted to complete their life cycle in saline soil keep the transcript expression of stress-responsive genes constitutively higher in the optimum growth environments, which can be further increased by several folds under stress conditions. The transcript expression of SbNHX1 gene, cloned from a leafless succulent halophyte Salicornia brachiata, was up-regulated under salinity stress, but its transcriptional regulation has not been studied so far. In the present study, a 1727\\u202fbp putative promoter (upstream to translation start site) of the SbNHX1 gene was cloned using a genome walking method. The bioinformatics analysis identified important stress-responsive cis-regulatory motifs, GT1, MBS, LTR and ARE, in addition to two leaf-specific enhancer motifs. The GUS expression analysis of stable transgenic tobacco plants, transformed with a transcriptional fusion of GUS with the full SbNHX1 promoter (NP1) or any of its five deletion fragments (NP2 to NP6), showed that the deletion of two enhancer motifs resulted in the sudden decrease in GUS expression in leaves but not in the stem or root tissues. In contrast, under salinity stress, the higher induction of GUS expression observed in NP1 and NP2 was correlated by the presence of salt-inducible GT1- and MBS-motifs which is distributed only in NP1 and NP2 deletion promoter fragments. Finally, we concluded that the SbNHX1 promoter has a 624\\u202fbp (-1727 to -1103\\u202fbp) regulatory region which contains the two leaf-specific enhancer motifs and salinity stress-inducible GT-1 and MBS motifs. We suggest the SbNHX1 gene promoter and fragments as a candidate alternative promoter/s for crop engineering for better stress tolerance, which can be amended according to the desired level of expression needed.", "sentence": [{"charStart": 0, "charEnd": 132, "index": 0}, {"index": 1, "charEnd": 283, "charStart": 134}, {"index": 2, "charEnd": 545, "charStart": 285}, {"index": 3, "charEnd": 756, "charStart": 547}, {"index": 4, "charEnd": 904, "charStart": 758}, {"index": 5, "charEnd": 1067, "charStart": 906}, {"index": 6, "charEnd": 1414, "charStart": 1069}, {"index": 7, "charEnd": 1652, "charStart": 1416}, {"index": 8, "charEnd": 1854, "charStart": 1654}, {"index": 9, "charEnd": 2066, "charStart": 1856}], "mesh": ["Chenopodiaceae", "Cloning, Molecular", "Gene Expression Regulation, Plant", "Plant Proteins", "Promoter Regions, Genetic", "Salinity", "Sequence Analysis, DNA", "Sodium-Hydrogen Exchangers", "Stress, Physiological", "Up-Regulation"], "type": ["Journal Article"]}', '{"docId": "30341356", "title": {"charStart": 0, "charEnd": 145}, "text": "Genotype\\u2009\\u00d7\\u2009Environment interactions of Nagina22 rice mutants for yield traits under low phosphorus, water limited and normal irrigated conditions. Multi environment testing helps identify stable genotypes especially for adverse abiotic stress situations. In the era of climate change and multiple abiotic stresses, it becomes important to analyze stability of rice lines under both irrigated and stress conditions. Mutants are an important genetic resource which can help in revealing the basis of natural variation. We analyzed 300 EMS induced mutants of aus rice cultivar Nagina22 (N22) for their G\\u2009\\u00d7\\u2009E interaction and stability under low phosphorus (P), water limited and irrigated conditions. Environmental effect and interaction were more significant than genotypic contribution on grain yield (GY), productive tillers (TN) and plant height (PH) under these three environmental conditions in dry season, 2010. GY and TN were more affected by low P stress than by water limited condition, but PH was not significantly different under these two stresses. Mutants G17, G209, G29, G91, G63 and G32 were stable for GY in decreasing order of stability across the three environments but G254 and G50 were stable only in low P, G17 and G45 only in water limited and G295 and G289 only in normal irrigated condition. We then selected and evaluated 3 high yielding mutants, 3 low yielding mutants and N22 for their stability and adaptability to these 3 environments in both wet and dry seasons for six years (2010-2015). The most stable lines based on the combined analysis of 12 seasons were G125 (NH210) under normal condition, G17 (NH686), G176 (NH363) and G284 (NH162) in low P condition and G176 (NH363) under water limited condition. G176 was the best considering all 3 conditions. When screened for 15 Pup1 gene-specific markers, G176 showed alleles similar to N22. While two other low-P tolerant lines G17 and G65 showed N22 similar alleles only at k-1 and k-5 but a different allele or null allele at 13 other loci. These stable mutants are a valuable resource for varietal development and to discover genes for tolerance to multiple abiotic stresses.", "sentence": [{"charStart": 0, "charEnd": 145, "index": 0}, {"index": 1, "charEnd": 253, "charStart": 147}, {"index": 2, "charEnd": 413, "charStart": 255}, {"index": 3, "charEnd": 515, "charStart": 415}, {"index": 4, "charEnd": 695, "charStart": 517}, {"index": 5, "charEnd": 913, "charStart": 697}, {"index": 6, "charEnd": 1056, "charStart": 915}, {"index": 7, "charEnd": 1311, "charStart": 1058}, {"index": 8, "charEnd": 1514, "charStart": 1313}, {"index": 9, "charEnd": 1733, "charStart": 1516}, {"index": 10, "charEnd": 1781, "charStart": 1735}, {"index": 11, "charEnd": 1866, "charStart": 1783}, {"index": 12, "charEnd": 2018, "charStart": 1868}, {"index": 13, "charEnd": 2154, "charStart": 2020}], "mesh": ["Agricultural Irrigation", "Edible Grain", "Environmental Exposure", "Gene-Environment Interaction", "Genotype", "Oryza", "Phosphorus"], "type": ["Journal Article", "Research Support, Non-U.S. Gov\'t"]}', '{"docId": "30342039", "title": {"charStart": 0, "charEnd": 100}, "text": "Combined effects of NaCl and Cd2+ stress on the photosynthetic apparatus of Thellungiella salsuginea. Plants show complex responses to abiotic stress while, the effect of the stress combinations can be different to those seen when each stress is applied individually. Here, we report on the effects of salt and/or cadmium on photosynthetic apparatus of Thellungiella salsuginea. Our results showed a considerable reduction of plant growth with some symptoms of toxicity, especially with cadmium treatment. The structural integrity of both photosystems (PSI and PSII) was mostly maintained under salt stress. Cadmium induced a considerable decrease of both PSI and PSII quantum yields and the electron transport rate ETR(I) and ETR(II) paralleled by an increase of non-photochemical quenching (NPQ). In addition, cadmium alone affects the rate of primary photochemistry by an increase of fluorescence at O-J phase and also the photo-electrochemical quenching at J-I phase. A positive L-band appeared with (Cd) treatment as an indicator of lower PSII connectivity, and a positive K-band reflecting the imbalance in number of electrons at donor and acceptor side. In continuity to our previous studies which showed that NaCl supply reduced Cd2+ uptake and limited its accumulation in shoot of divers halophyte species, here as a consequence, we demonstrated the NaCl-induced enhancement effect of Cd2+ toxicity on the PSII activity by maintaining the photosynthetic electron transport chain as evidenced by the differences in \\u03c8O, \\u03c6Eo, ABS/RC and TR0/RC and by improvement of performance index PI(ABS), especially after short time of treatment. A significant decrease of LHCII, D1 and CP47 amounts was detected under (Cd) treatment. However, NaCl supply alleviates the Cd2+ effect on protein abundance including LHCII and PSII core complex (D1 and CP47).", "sentence": [{"charStart": 0, "charEnd": 100, "index": 0}, {"index": 1, "charEnd": 266, "charStart": 102}, {"index": 2, "charEnd": 377, "charStart": 268}, {"index": 3, "charEnd": 504, "charStart": 379}, {"index": 4, "charEnd": 606, "charStart": 506}, {"index": 5, "charEnd": 797, "charStart": 608}, {"index": 6, "charEnd": 970, "charStart": 799}, {"index": 7, "charEnd": 1159, "charStart": 972}, {"index": 8, "charEnd": 1639, "charStart": 1161}, {"index": 9, "charEnd": 1727, "charStart": 1641}, {"index": 10, "charEnd": 1849, "charStart": 1729}], "mesh": ["Biomass", "Brassicaceae", "Cadmium", "Chlorophyll", "Electron Transport", "Fluorescence", "Photosynthesis", "Photosystem I Protein Complex", "Photosystem II Protein Complex", "Sodium Chloride", "Stress, Physiological", "Thylakoids"], "type": ["Journal Article", "Research Support, Non-U.S. Gov\'t"]}', '{"docId": "30342485", "title": {"charStart": 0, "charEnd": 109}, "text": "Shared and genetically distinct Zea mays transcriptome responses to ongoing and past low temperature exposure. BACKGROUND: Cold temperatures and their alleviation affect many plant traits including the abundance of protein coding gene transcripts. Transcript level changes that occur in response to cold temperatures and their alleviation are shared or vary across genotypes. In this study we identify individual transcripts and groups of functionally related transcripts that consistently respond to cold and its alleviation. Genes that respond differently to temperature changes across genotypes may have limited functional importance. We investigate if these genes share functions, and if their genotype-specific gene expression levels change in magnitude or rank across temperatures. RESULTS: We estimate transcript abundances from over 22,000 genes in two unrelated Zea mays inbred lines during and after cold temperature exposure. Genotype and temperature contribute to many genes\' abundances. Past cold exposure affects many fewer genes. Genes up-regulated in cold encode many cytokinin glucoside biosynthesis enzymes, transcription factors, signalling molecules, and proteins involved in diverse environmental responses. After cold exposure, protease inhibitors and cuticular wax genes are newly up-regulated, and environmentally responsive genes continue to be up-regulated. Genes down-regulated in response to cold include many photosynthesis, translation, and DNA replication associated genes. After cold exposure, DNA replication and translation genes are still preferentially downregulated. Lignin and suberin biosynthesis are newly down-regulated. DNA replication, reactive oxygen species response, and anthocyanin biosynthesis genes have strong, genotype-specific temperature responses. The ranks of genotypes\' transcript abundances often change across temperatures. CONCLUSIONS: We report a large, core transcriptome response to cold and the alleviation of cold. In cold, many of the core suite of genes are up or downregulated to control plant growth and photosynthesis and limit cellular damage. In recovery, core responses are in part to prepare for future stress. Functionally related genes are consistently and greatly up-regulated in a single genotype in response to cold or its alleviation, suggesting positive selection has driven genotype-specific temperature responses in maize.", "sentence": [{"charStart": 0, "charEnd": 109, "index": 0}, {"index": 1, "charEnd": 246, "charStart": 111}, {"index": 2, "charEnd": 374, "charStart": 248}, {"index": 3, "charEnd": 525, "charStart": 376}, {"index": 4, "charEnd": 636, "charStart": 527}, {"index": 5, "charEnd": 786, "charStart": 638}, {"index": 6, "charEnd": 935, "charStart": 788}, {"index": 7, "charEnd": 998, "charStart": 937}, {"index": 8, "charEnd": 1043, "charStart": 1000}, {"index": 9, "charEnd": 1227, "charStart": 1045}, {"index": 10, "charEnd": 1382, "charStart": 1229}, {"index": 11, "charEnd": 1503, "charStart": 1384}, {"index": 12, "charEnd": 1602, "charStart": 1505}, {"index": 13, "charEnd": 1660, "charStart": 1604}, {"index": 14, "charEnd": 1800, "charStart": 1662}, {"index": 15, "charEnd": 1880, "charStart": 1802}, {"index": 16, "charEnd": 1977, "charStart": 1882}, {"index": 17, "charEnd": 2112, "charStart": 1979}, {"index": 18, "charEnd": 2182, "charStart": 2114}, {"index": 19, "charEnd": 2403, "charStart": 2184}], "mesh": ["Cold Temperature", "Environment", "Gene Expression Profiling", "Genotype", "Glucose", "Photosynthesis", "RNA, Messenger", "Signal Transduction", "Transcription, Genetic", "Up-Regulation", "Zea mays"], "type": ["Journal Article"]}', '{"docId": "30347736", "title": {"charStart": 0, "charEnd": 128}, "text": "Genome-Wide Characterization of the sHsp Gene Family in Salix suchowensis Reveals Its Functions under Different Abiotic Stresses. Small heat shock proteins (sHsps) function mainly as molecular chaperones that play vital roles in response to diverse stresses, especially high temperature. However, little is known about the molecular characteristics and evolutionary history of the sHsp family in Salix suchowensis, an important bioenergy woody plant. In this study, 35 non-redundant sHsp genes were identified in S. suchowensis, and they were divided into four subfamilies (C, CP, PX, and MT) based on their phylogenetic relationships and predicted subcellular localization. Though the gene structure and conserved motif were relatively conserved, the sequences of the Hsp20 domain were diversified. Eight paralogous pairs were identified in the Ssu-sHsp family, in which five pairs were generated by tandem duplication events. Ka/Ks analysis indicated that Ssu-sHsps had undergone purifying selection. The expression profiles analysis showed Ssu-Hsps tissue-specific expression patterns, and they were induced by at least one abiotic stress. The expression correlation between two paralogous pairs (Ssu-sHsp22.2-CV/23.0-CV and 23.8-MT/25.6-MT) were less than 0.6, indicating that they were divergent during the evolution. Various cis-acting elements related to stress responses, hormone or development, were detected in the promoter of Ssu-sHsps. Furthermore, the co-expression network revealed the potential mechanism of Ssu-sHsps under stress tolerance and development. These results provide a foundation for further functional research on the Ssu-sHsp gene family in S. suchowensis.", "sentence": [{"charStart": 0, "charEnd": 128, "index": 0}, {"index": 1, "charEnd": 286, "charStart": 130}, {"index": 2, "charEnd": 449, "charStart": 288}, {"index": 3, "charEnd": 673, "charStart": 451}, {"index": 4, "charEnd": 798, "charStart": 675}, {"index": 5, "charEnd": 926, "charStart": 800}, {"index": 6, "charEnd": 1001, "charStart": 928}, {"index": 7, "charEnd": 1141, "charStart": 1003}, {"index": 8, "charEnd": 1321, "charStart": 1143}, {"index": 9, "charEnd": 1446, "charStart": 1323}, {"index": 10, "charEnd": 1571, "charStart": 1448}, {"index": 11, "charEnd": 1685, "charStart": 1573}], "mesh": ["Evolution, Molecular", "Heat-Shock Proteins", "Heat-Shock Response", "Phylogeny", "Plant Proteins", "Salix", "Salt Stress"], "type": ["Journal Article"]}', '{"docId": "30348322", "title": {"charStart": 0, "charEnd": 105}, "text": "Overexpression of Arabidopsis ubiquitin ligase AtPUB46 enhances tolerance to drought and oxidative stress. The U-Box E3 ubiquitin ligase, AtPUB46, functions in the drought response: T-DNA insertion mutants of this single paralogous gene are hypersensitive to water- and oxidative stress (Adler et al. BMC Plant Biology 17:8, 2017). Here we analyze the phenotype of AtPUB46 overexpressing (OE) plants. AtPUB46-OE show increased tolerance to water stress and have smaller leaf blades and reduced stomatal pore area and stomatal index compared with wild type (WT). Despite this, the rate of water loss from detached rosettes is similar in AtPUB46-OE and WT plants. Germination of AtPUB46-OE seeds was less sensitive to salt than WT whereas seedling greening was more sensitive. We observed a complex response to oxidative stress applied by different agents: AtPUB46-OE plants were hypersensitive to H2O2 but hyposensitive to methyl viologen. AtPUB46-GFP fusion protein is cytoplasmic, however, in response to H2O2 a considerable proportion translocates to the nucleus. We conclude that the differential stress phenotype of the AtPUB46-OE does not result from its smaller leaf size but from a change in the activity of a stress pathway(s) regulated by a degradation substrate of the AtPUB46 E3 and also from a reduction in stomatal pore size and index.", "sentence": [{"charStart": 0, "charEnd": 105, "index": 0}, {"index": 1, "charEnd": 299, "charStart": 107}, {"index": 2, "charEnd": 330, "charStart": 301}, {"index": 3, "charEnd": 399, "charStart": 332}, {"index": 4, "charEnd": 560, "charStart": 401}, {"index": 5, "charEnd": 660, "charStart": 562}, {"index": 6, "charEnd": 773, "charStart": 662}, {"index": 7, "charEnd": 937, "charStart": 775}, {"index": 8, "charEnd": 1064, "charStart": 939}, {"index": 9, "charEnd": 1347, "charStart": 1066}], "mesh": ["Arabidopsis", "Arabidopsis Proteins", "Cytoplasm", "Dehydration", "Droughts", "Genes, Reporter", "Germination", "Hydrogen Peroxide", "Oxidative Stress", "Plant Leaves", "Plant Stomata", "Recombinant Fusion Proteins", "Seedlings", "Seeds", "Sodium Chloride", "Stress, Physiological", "Ubiquitin-Protein Ligases"], "type": ["Journal Article"]}', '{"docId": "30348330", "title": {"charStart": 0, "charEnd": 89}, "text": "Low-temperature tolerance in land plants: Are transcript and membrane responses conserved? Plants\' tolerance of low temperatures is an economically and ecologically important limitation on geographic distributions and growing seasons. Tolerance for low temperatures varies significantly across different plant species, and different mechanisms likely act in different species. In order to survive low-temperature stress, plant membranes must maintain their fluidity in increasingly cold and oxidative cellular environments. The responses of different species to low-temperature stress include changes to the types and desaturation levels of membrane lipids, though the precise lipids affected tend to vary by species. Regulation of membrane dynamics and other low-temperature tolerance factors are controlled by both transcriptional and post-transcriptional mechanisms. Here, we review low-temperature induced changes in both membrane lipid composition and gene transcription across multiple related plant species with differing degrees of low-temperature tolerance. We attempt to define a core set of changes for transcripts and lipids across species and treatment variations. Some responses appear to be consistent across all species for which data are available, while many others appear likely to be species or family-specific. Potential rationales are presented, including variance in testing, reporting and the importance of considering the level of stress perceived by the plant.", "sentence": [{"charStart": 0, "charEnd": 89, "index": 0}, {"index": 1, "charEnd": 233, "charStart": 91}, {"index": 2, "charEnd": 375, "charStart": 235}, {"index": 3, "charEnd": 522, "charStart": 377}, {"index": 4, "charEnd": 716, "charStart": 524}, {"index": 5, "charEnd": 868, "charStart": 718}, {"index": 6, "charEnd": 1065, "charStart": 870}, {"index": 7, "charEnd": 1176, "charStart": 1067}, {"index": 8, "charEnd": 1330, "charStart": 1178}, {"index": 9, "charEnd": 1485, "charStart": 1332}], "mesh": ["Acclimatization", "Cold Temperature", "Embryophyta", "Membrane Lipids", "Species Specificity", "Stress, Physiological", "Transcriptome"], "type": ["Journal Article", "Review"]}', '{"docId": "30350236", "title": {"charStart": 0, "charEnd": 150}, "text": "Genome-wide analysis of genes encoding MBD domain-containing proteins from tomato suggest their role in fruit development and abiotic stress responses. In tomato, DNA methylation has an inhibitory effect on fruit ripening. The inhibition of DNA methyltransferase by 5-azacytidine results in premature fruit ripening. Methyl CpG binding domain (MBD) proteins are the readers of DNA methylation marks and help in the recruitment of chromatin-modifying enzymes which affect gene expression. Therefore, we investigate their contribution during fruit development. In this study, we identified and analyzed 18 putative genes of Solanum lycopersicum and Solanum pimpinellifolium encoding MBD proteins. We also identified tomato MBD syntelogs in Capsicum annum and Solanum tuberosum. Sixty-three MBD genes identified from four different species of solanaceae were classified into three groups. An analysis of the conserved domains in these proteins identified additional domains along with MBD motif. The transcript profiling of tomato MBDs in wild-type and two non-ripening mutants, rin and Nr, indicated constructive information regarding their involvement during fruit development. When we performed a stage-specific expression analysis during fruit ripening, a gradual decrease in transcript accumulation in the wild-type fruit was detected. However, a very low expression was observed in the ripening mutants. Furthermore, many ethylene-responsive cis-elements were found in SlMBD gene promoters, and some of them were induced in the presence of exogenous ethylene. Further, we detected the possible role of these MBDs in abiotic stresses. We found that few genes were differentially expressed under various abiotic stress conditions. Our results provide an evidence of the involvement of the tomato MBDs in fruit ripening and abiotic stress responses, which would be helpful in further studies on these genes in tomato fruit ripening.", "sentence": [{"charStart": 0, "charEnd": 150, "index": 0}, {"index": 1, "charEnd": 221, "charStart": 152}, {"index": 2, "charEnd": 315, "charStart": 223}, {"index": 3, "charEnd": 486, "charStart": 317}, {"index": 4, "charEnd": 557, "charStart": 488}, {"index": 5, "charEnd": 693, "charStart": 559}, {"index": 6, "charEnd": 774, "charStart": 695}, {"index": 7, "charEnd": 884, "charStart": 776}, {"index": 8, "charEnd": 991, "charStart": 886}, {"index": 9, "charEnd": 1175, "charStart": 993}, {"index": 10, "charEnd": 1336, "charStart": 1177}, {"index": 11, "charEnd": 1405, "charStart": 1338}, {"index": 12, "charEnd": 1561, "charStart": 1407}, {"index": 13, "charEnd": 1635, "charStart": 1563}, {"index": 14, "charEnd": 1730, "charStart": 1637}, {"index": 15, "charEnd": 1931, "charStart": 1732}], "mesh": ["Capsicum", "DNA Methylation", "Ethylenes", "Fruit", "Gene Expression Profiling", "Gene Expression Regulation, Plant", "Genes, Plant", "Genome-Wide Association Study", "Lycopersicon esculentum", "Methyl CpG Binding Domain", "Plant Proteins", "Promoter Regions, Genetic", "Solanum tuberosum", "Stress, Physiological"], "type": ["Journal Article"]}']
    extract(parts)