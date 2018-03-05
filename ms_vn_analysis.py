# Script for exporting audubon core CSV file representing UF specimens.
#
# Author: Julie Winchester <julia.m.winchester@gmail.com>
# February 14, 2018

from requests import Request, Session
import credentials
import csv
import pymysql
import requests
import cPickle as pickle

def db_conn():
	return pymysql.connect(unix_socket = credentials.db['socket'],
						   user = credentials.db['username'],
						   password = credentials.db['password'],
						   db = credentials.db['db'],
						   charset = 'utf8mb4',
						   cursorclass=pymysql.cursors.DictCursor)

def db_query(cursor, sql, args=None):
	if args is not None:
		args = [args]
	cursor.execute(sql, args)
	return cursor.fetchall()

def vn_url(genus, species):
	return 'http://api.vertnet-portal.appspot.com/api/search?q=%7B%22q%22:%22genus:' + genus + '%20specificepithet:' + species + '%22%7D'

def create_tstr(kingdom='', phylum='', classs='', order='', family='', genus=''):
	a = [kingdom.capitalize(), phylum.capitalize(), classs.capitalize(), order.capitalize(), family.capitalize(), genus.capitalize()]
	b = [x.replace('/', '') for x in a]
	c = [x.replace(' ', '') for x in b]
	return '/'.join(c)

def split_tstr(tstr):
	return tstr.split('/')

def validate_tstr(tstr):
	if tstr.count('/') == 5:
		return True
	else:
		return False

def cmp_str_loose(a, b):
	if a == b:
		return 1
	elif a == '' or b == '':
		return 2
	else:
		return 0

def mask_loose(a, b):
	x = []
	for i in range(0, len(a)):
		x.append(cmp_str_loose(a[i], b[i]))
	return x

def compare_partial_str(tstr1, tstr2):
	mask = mask_loose(tstr1, tstr2)
	if 0 in mask:
		return 0
	elif 2 in mask:
		return 2
	else:
		return 1

def compare_tstr(tstr1, tstr2):
	if not validate_tstr(tstr1) or not validate_tstr(tstr2):
		raise ValueError('Problem validating tstr ' +str(tstr1) + ' or ' + str(tstr2))
	if tstr1 == tstr2:
		return 1
	else:
		tstr1_array = split_tstr(tstr1)
		tstr2_array = split_tstr(tstr2)
		comp = compare_partial_str(tstr1_array, tstr2_array)
		if comp == 1:
			return 1
		elif comp == 2:
			return 2
		else:
			return 0

def merge_tstr(tstr1, tstr2):
	tstr1_array = split_tstr(tstr1)
	tstr2_array = split_tstr(tstr2)
	mask = mask_loose(tstr1_array, tstr2_array)
	if 0 in mask or 2 not in mask:
		return False
	else:
		x = []
		for i in range(0,len(tstr1_array)):
			if mask[i] == 1:
				x.append(tstr1_array[i])
			elif mask[i] == 2:
				if tstr1_array[i] == '' and tstr2_array[i] != '':
					x.append(tstr2_array[i])
				elif tstr2_array[i] == '' and tstr1_array[i] != '':
					x.append(tstr1_array[i])
				else:
					x.append('')
		return x

# Get MorphoSource taxonomies

conn = db_conn()
c = conn.cursor()

sql = """ SELECT DISTINCT genus, species FROM ms_taxonomy_names """

r = db_query(c, sql)

# Create genus dict

genusDict = {}
for t in r:
	pickle.dump(genusDict, open("genusDict.p", "wb"))
	print t['genus'] + ' ' + t['species']
	if t['genus'] not in genusDict:
		genusDict[t['genus']] = {}
		genusDict[t['genus']]['sp'] = []
		genusDict[t['genus']]['taxonomies'] = []
	genusDict[t['genus']]['sp'].append(t['species'])
	resp = requests.get(vn_url(t['genus'], t['species']))
	try:
		json = resp.json()
		for specimen in json['recs']:
			if 'kingdom' in specimen:
				k = specimen['kingdom']
			else:
				k = ''
			if 'phylum' in specimen:
				p = specimen['phylum']
			else:
				p = ''
			if 'class' in specimen:
				c = specimen['class']
			else:
				c = ''
			if 'order' in specimen:
				o = specimen['order']
			else:
				o = ''
			if 'family' in specimen:
				f = specimen['family']
			else:
				f = ''
			if 'genus' in specimen:
				g = specimen['genus']
			else:
				g = ''
			if 'specificepithet' in specimen:
				s = specimen['specificepithet']
			else:
				s = ''
			genusDict[t['genus']]['taxonomies'].append({
				'kingdom' : k,
				'phylum' : p,
				'class' : c,
				'order' : o,
				'family' : f,
				'genus' : g,
				'species' : s 
			})
	except:
		continue	

pickle.dump(genusDict, open("genusDict.p", "wb"))

# Create unique tstr dict
for genus, d in genusDict.iteritems():
	genusDict[genus]['uniq_tn'] = {}
	for tn in genusDict[genus]['taxonomies']:
		tstr_new = create_tstr(tn['kingdom'], tn['phylum'], tn['class'], tn['order'], tn['family'], tn['genus'])
		if len(genusDict[genus]['uniq_tn'].keys()) == 0:
			genusDict[genus]['uniq_tn'][tstr_new] = 1
		else:
			is_new = 1
			to_be_merged = 0
			merge_keys = []
			for tstr in genusDict[genus]['uniq_tn'].keys():
				comp = compare_tstr(tstr, tstr_new)
				if comp == 1:
					is_new = 0
					genusDict[genus]['uniq_tn'][tstr] += 1
					break
				elif comp == 2:
					is_new = 0
					to_be_merged = 1
					merge_keys.append(tstr)
					break 
			if is_new:
				genusDict[genus]['uniq_tn'][tstr_new] = 1
			if to_be_merged:
				del_array = []
				for k in merge_keys:
					m = merge_tstr(k, tstr_new)
					new_key = create_tstr(*m)
					genusDict[genus]['uniq_tn'][new_key] = genusDict[genus]['uniq_tn'][k] + 1
					del_array.append(k)
				for k in del_array:
					del genusDict[genus]['uniq_tn'][k]

pickle.dump(genusDict, open("genusDict.p", "wb"))

# Write output CSVs

full_out = open('all_unique_vn_taxonomies.csv', 'w')
with full_out:
	writer = csv.writer(full_out)
	for genus, d in genusDict.iteritems():
		writer.writerow([genus])
		for tn in genusDict[genus]['uniq_tn'].keys():
			tn_array = tn.split('/')
			tn_array.append(genusDict[genus]['uniq_tn'][tn])
			writer.writerow(tn_array)
		writer.writerow([''])

sum_out = open('unique_vn_taxonomies_count.csv', 'w')
with sum_out:
	writer = csv.writer(sum_out)
	writer.writerow(['genus', 'number of taxonomies', 'proportions'])
	for genus, d in genusDict.iteritems():
		num_tn = len(d['uniq_tn'].keys())
		tn_counts = d['uniq_tn'].values()
		sort_tn_counts = sorted(tn_counts, reverse=True)
		print sort_tn_counts
		tn_counts_sum = sum(sort_tn_counts)
		sort_tn_prop = [float(x)/float(tn_counts_sum) for x in sort_tn_counts]
		print sort_tn_prop
		genus_array = sort_tn_prop
		genus_array.insert(0, num_tn)
		genus_array.insert(0, genus)
		writer.writerow(genus_array)





