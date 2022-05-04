from pmid_extractor import extract
json_file = "pmids.json"


def load_json():
    """
    This function loads the input json file and returns the
    json data
    :param file_name: the JSON file that needs to be loaded
    :return: the JSON data
    """
    json_file = input("Please enter the JSON file name: ")
    if json_file is not None and json_file.endswith(".json"):
        with open(json_file) as json_obj:
            # get the json data from the file
            json_file_content = json_obj.read().split('\n')
            print(f"There are {len(json_file_content)} PMIDs.")
            for i in range(0, len(json_file_content), 20):
                # split the list into parts of 20 pmids.
                parts = json_file_content[i:i+20]
                # now, add these into the main script
                extract(parts)
    else:
        print("Invalid JSON file.")


if __name__ == '__main__':
    load_json()