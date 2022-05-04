import pymongo

mongo_client = pymongo.MongoClient("mongodb://localhost:27017/")
db = mongo_client['database']
collection = db['pmid_output']
# test records
records = [{"name": "user1", "hobby": "painting"}, {"name": "user2", "hobby": "singing"}]


def save(records):
    collection.insert_many(records)
    print("Successfully inserted the records into database")
