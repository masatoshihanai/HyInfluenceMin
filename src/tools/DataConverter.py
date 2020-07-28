import datetime
import pickle
import sys
import argparse
from functools import cmp_to_key

newCategories = [
    'Food & Restaurant',
    'School',
    'Transportation',
    'Office',
    'Shops & Stores',
    'Outdoor Misc',
    'Indoor Misc',
    'Others'
    ]
newCategoriesIndexMap = {}
for i in range(len(newCategories)):
    newCategoriesIndexMap[newCategories[i]] = i
#print(newCategoriesIndexMap)

categoriesDict = {}
f = open('categories-map.txt')
while True:
    line = f.readline()
    if not line:
        break
    line = line.strip()
    if line[0] == '#':
        continue
    pair = line.split(',')
    if len(pair) != 2:
        print('Mapping error: {}'.format(pair))
    else:
        key, value = pair
        categoriesDict[key] = value
f.close()

def utc2iso(str):
    return datetime.datetime.strptime(str, '%a %b %d %H:%M:%S %z %Y').strftime('%Y-%m-%dT%H:%M:%SZ')

class Record:
    def __init__(self, datapack):
        self.userId = int(datapack[0])
        self.venueId = int(datapack[1], 16)
        self.venueCategoryId = int(datapack[2], 16)
        self.venueCategory = datapack[3]
        self.convCategory = categoriesDict[self.venueCategory]
        self.convCategoryId = newCategoriesIndexMap[self.convCategory]
        self.latitude = float(datapack[4])
        self.longitude = float(datapack[5])
        self.timezoneOffset = datapack[6]
        self.utcTimestamp = datapack[7]
        self.isoTimestamp = utc2iso(datapack[7])
    def comparator(a, b):
        if a.isoTimestamp < b.isoTimestamp:
            return -1
        elif a.isoTimestamp > b.isoTimestamp:
            return 1
        else:
            return 0
    def __str__(self):
        return '{{userId: {}, venueId: {}, venueCategoryId: {}, venueCategory: {}, latitude: {}, longitude: {}, timezoneOffset: {}, utcTimestamp: {}, isoTimestamp: {}}}'.format(
            self.userId, self.venueId, self.venueCategoryId, self.venueCategory, self.latitude, self.longitude, self.timezoneOffset, self.utcTimestamp, self.isoTimestamp)

def loadData(path):
    print('Reading {}'.format(path))
    f = open(path)

    data_complete = {}
    data_complete['checkin'] = 0
    data_complete['users'] = 0
    data_complete['locations'] = 0
    data_complete['categories'] = 0
    data_complete['conv_categories'] = 0
    data_complete['usr_map'] = {}
    data_complete['loc_map'] = {}
    data_complete['cat_map'] = {}
    data_complete['data'] = {}

    data = {}   # usrID: dataObj

    usr_ct = 0
    loc_ct = 0
    cat_ct = 0

    categorySet = set()
    convCategorySet = set()
    locations = set()
    categoryIDs = set()

    while True:
        line = f.readline()
        if not line:
            break
        datapack = line.strip().split(',')
        if datapack[0] == 'userId':
            continue

        #print(datapack)
        dataObj = Record(datapack)
        categorySet.add(dataObj.venueCategory)
        convCategorySet.add(dataObj.convCategory)
        #print(dataObj)

        if dataObj.userId not in data:
            data[dataObj.userId] = []
            data_complete['users'] += 1
            data_complete['usr_map'][dataObj.userId] = usr_ct
            usr_ct += 1
        data[dataObj.userId].append(dataObj)

        if dataObj.venueId not in locations:
            locations.add(dataObj.venueId)
            data_complete['loc_map'][dataObj.venueId] = loc_ct
            loc_ct += 1

        if dataObj.venueCategoryId not in categoryIDs:
            categoryIDs.add(dataObj.venueCategoryId)
            data_complete['cat_map'][dataObj.venueCategoryId] = cat_ct
            cat_ct += 1

        data_complete['checkin'] += 1
        sys.stderr.write('\rLoaded: {}, Categories: {}, Users: {}, Locations: {}'.format(data_complete['checkin'], len(categorySet), data_complete['users'], len(locations)))
        sys.stderr.flush()

    f.close()

    print('')
    print('Sorting Category Set...')
    categorySet = sorted(categorySet)

    print('Normalizing Data...')
    for usrId in data:
        mapped_usrId = data_complete['usr_map'][usrId]
        for j in range(len(data[usrId])):
            data[usrId][j].userId = mapped_usrId
            data[usrId][j].venueId = data_complete['loc_map'][data[usrId][j].venueId]
            data[usrId][j].venueCategoryId = data_complete['cat_map'][data[usrId][j].venueCategoryId]
        data_complete['data'][mapped_usrId] = data[usrId]

    print('Sorting Data...')
    for usrId in data_complete['data']:
        data_complete['data'][usrId] = sorted(data_complete['data'][usrId], key=cmp_to_key(Record.comparator))
    data_complete['data'] = sorted(data_complete['data'].items())

    data_complete['locations'] = len(locations)
    data_complete['categories'] = len(categorySet)
    data_complete['conv_categories'] = len(convCategorySet)

    return data_complete, categorySet

def dump(data, path):
    pickle.dump(data, open(path, 'wb'))
    print('data dumped as {}'.format(path))

def load(path):
    data = pickle.load(open(path, 'rb'))
    print('data loaded from {}'.format(path))
    return data

def writeData(data, categorySet, path):
    cat_path = '{}_categories.txt'.format(path)
    print('Recording categories to {}'.format(cat_path))
    f = open(cat_path, 'w')
    f.write('# Categories: {}\n'.format(len(categorySet)))
    for category in categorySet:
        f.write('{}\n'.format(category))
    f.close()

    conv_path = '{}.checkin'.format(path)
    f = open(conv_path, 'w')

    f.write('## # CheckIn: {}. # Users: {}. # Locations: {}. # Original Categories: {}. # Categories: {}\n'.format(data['checkin'], data['users'], data['locations'],
                                                                                                                   data['categories'], data['conv_categories']))
    f.write('## [user]\t[check-in time]\t[latitude]\t[longitude]\t[location id]\t[orig_category]\t[orig_category id]\t[category]\t[category id]\n')
    records = data['data']
    counter = 0
    for pair in records:
        record_list = pair[1]
        for record in record_list:
            f.write('{},{},{},{},{},{},{},{},{}\n'.format(record.userId, record.isoTimestamp, record.latitude, record.longitude, record.venueId,
                                                          record.venueCategory, record.venueCategoryId,
                                                          record.convCategory, record.convCategoryId))
            counter += 1
            sys.stderr.write('\rConverting data to {}: {}/{}'.format(conv_path, counter, data['checkin']))
            sys.stderr.flush()

    f.close()
    print('\nConvertion Finished.')

if __name__ == '__main__':
    ### Usage example
    # py DataConverter.py -d data/org/dataset_TSMC2014_NYC.csv -s dataset_TSMC2014_NYC -o dataset_TSMC2014_NYC
    # py DataConverter.py -d data/org/dataset_TSMC2014_TKY.csv -s dataset_TSMC2014_TKY -o dataset_TSMC2014_TKY
    # py DataConverter.py -l dataset_TSMC2014_NYC -o dataset_TSMC2014_NYC
    # py DataConverter.py -l dataset_TSMC2014_TKY -o dataset_TSMC2014_TKY

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--save', type=str, default=None,
                        help='save/dump to cwd? If save, specify a prefix name; If not, None')
    parser.add_argument('-l', '--load', type=str, default=None,
                        help='load from cwd? If load, specify a prefix name; If not, None')
    parser.add_argument('-d', '--dataset', type=str, default=None,
                        help='If not load, specify the path of the dataset')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='prefix name for the output results (will output to the cwd)')
    args = parser.parse_args()

    data, categorySet = None, None

    if not args.output:
        print('Need to specify the output prefix name! Try "-h" for more details')
        exit(-1)

    if args.dataset:
        data, categorySet = loadData(args.dataset)
    elif args.load:
        data = load('{}_obj.pk'.format(args.load))
        categorySet = load('{}_categories.pk'.format(args.load))
    else:
        print('Need to specify a dataset or a dumped data loading position! Try "-h" for more details')
        exit(-2)

    if args.save:
        dump(data, '{}_obj.pk'.format(args.save))
        dump(categorySet, '{}_categories.pk'.format(args.save))

    writeData(data, categorySet, args.output)