# Script for exporting audubon core CSV file representing UF specimens.
#
# Author: Julie Winchester <julia.m.winchester@gmail.com>
# February 14, 2018

from functools import reduce
from requests import Request, Session
import credentials
import csv
import operator
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

# Taxonomy string (tstr) functions

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

# Tree functions

def search_tree(t, value, rank=None, prepath=[], results=[]):
    for k, v in t.iteritems():
        path = prepath + [k]
        if k == value and (rank is None or rank == len(path)): # found value
            results.append(path)
            return results
        elif hasattr(v, 'items'): # v is a dict
            results = search_tree(v, value, rank, path, results) # recursive call
    return copy.deepcopy(results)

def get_from_dict(d, key_list):
    return reduce(operator.getitem, key_list, d)

def set_in_dict(d, key_list, value):
    get_from_dict(d, key_list[:-1])[key_list[-1]] = value


def merge_tree_recursive(t, src, tgt, path_end=[], del_paths=[]):
	del_paths = copy.deepcopy(del_paths)
	src_dict = get_from_dict(t, src+path_end)
	tgt_dict = get_from_dict(t, tgt+path_end)
	for k,v in src_dict.iteritems():
		if k in tgt_dict:
			if type(v) is dict:
				path_end.append(k)
				t, path_end = merge_tree_recursive(t, src, tgt, path_end, del_paths)
				path_end = path_end[:-1]
			elif type(v) is list:
				new_v = src_dict[k] + tgt_dict[k]
				tgt_dict[k] = new_v
			else:
				raise ValueError('Non-list terminal node at ' + str(src+path_end+[k]))
		else:
			tgt_dict[k] = v
	return t, path_end

def merge_tree_path(t, src, tgt):
	# Takes in a tree and two different paths, merges contents of src to target
	if src == tgt or len(src) != len(tgt) or src[-1] != tgt[-1]:
		return False
	t, path_end = merge_tree_recursive(t, src, tgt, [])
	del get_from_dict(t, src[:-1])[src[-1]]
	return t

# Misc functions

def get_tn_uniques(t):
	result = {'k': [], 'p': [], 'c': [], 'o': [], 'f': [], 'g': []}
	for k, k_val in t.iteritems():
		if k != '' and k not in result['k']:
			result['k'].append(k)
		for p, p_val in k_val.iteritems():
			if p != '' and p not in result['p']:
				result['p'].append(p)
			for c, c_val in p_val.iteritems():
				if c != '' and c not in result['c']:
					result['c'].append(c)
				for o, o_val in c_val.iteritems():
					if o != '' and o not in result['o']:
						result['o'].append(o)
					for f, f_val in o_val.iteritems():
						if f != '' and f not in result['f']:
							result['f'].append(f)
						for g, g_val in f_val.iteritems():
							if g != '' and g not in result['g']:
								result['g'].append(g)
	return result

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

# Create canonical taxonomy
for genus, d in genusDict.iteritems():
	if len(d['uniq_tn'].keys()) == 0:
		d['canon_tn'] = None
	else: 
		d['canon_tn'] = sorted(d['uniq_tn'].items(), key=operator.itemgetter(1), reverse=True)[0][0]

# Create taxonomy tree
t = {}
unsorted_genera = []
for genus, d in genusDict.iteritems():
	if not d['canon_tn']:
		unsorted_genera.append(genus)
	else:
		k, p, c, o, f, g = split_tstr(d['canon_tn'])
		if k not in t:
			t[k] = {}
		if p not in t[k]:
			t[k][p] = {}
		if c not in t[k][p]:
			t[k][p][c] = {}
		if o not in t[k][p][c]:
			t[k][p][c][o] = {}
		if f not in t[k][p][c][o]:
			t[k][p][c][o][f] = {}
		if g not in t[k][p][c][o][f]:
			t[k][p][c][o][f][g] = d['sp']

# output taxonomy tree
tree_out = open('canon_taxonomy_tree.csv', 'w')
with tree_out:
	writer = csv.writer(tree_out)
	writer.writerow(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species list'])
	for k in t.keys():
		writer.writerow([k])
		for p in t[k].keys():
			writer.writerow(['', p])
			for c in t[k][p].keys():
				writer.writerow(['', '', c])
				for o in t[k][p][c].keys():
					writer.writerow(['', '', '', o])
					for f in t[k][p][c][o].keys():
						writer.writerow(['', '', '', '', f])
						for g in t[k][p][c][o][f].keys():
							genus_array = ['', '', '', '', '', g]
							if g in genusDict:
								genus_array.extend(genusDict[g]['sp'])
							writer.writerow(genus_array)



# Associate hanging taxonomies with (any) fully fleshed out taxonomy
for k in t.keys():
	if k == '':
		# No kingdom info
		for p in t[k].keys():
			if p == '':
				# No kingdom, phylum info
				for c in t[k][p].keys():
					if c == '':
						# No kingdom, phylum, class info
						for o in t[k][p][c].keys():
							if o == '':
								# No kingdom, phylum, class, order info
								for f in t[k][p][c][o].keys():
									# Hanging family
									matches = search_tree(t, f, None, [], [])
									matches = [x for x in matches if x != [k, p, c, o, f]]
									print matches
									for i in range(0,len(matches)):
										if i == 0:
											t = merge_tree_path(t, [k, p, c, o, f], matches[i])
										else:
											t = merge_tree_path(t, matches[i-1], matches[i])

# For each family

tn_uniq = get_tn_uniques(t)

# Try and merge taxonomies based on sub-category voting (# of sub-categories)
for f in tn_uniq['f']:
	matches = search_tree(t, f, 5, [], [])
	match_len = {}
	for i in range(0, len(matches)):
		match_len[i] = len(get_from_dict(t, matches[i]))
	max_match_i = sorted(match_len.items(), key=operator.itemgetter(1), reverse=True)[0][0]
	match_max = matches[max_match_i]
	matches_minus_max = [x for x in matches if x != match_max]
	for m in matches_minus_max:
		t = merge_tree_path(t, m, match_max)

for o in tn_uniq['o']:
	matches = search_tree(t, o, 4, [], [])
	match_len = {}
	for i in range(0, len(matches)):
		match_len[i] = len(get_from_dict(t, matches[i]))
	max_match_i = sorted(match_len.items(), key=operator.itemgetter(1), reverse=True)[0][0]
	match_max = matches[max_match_i]
	matches_minus_max = [x for x in matches if x != match_max]
	for m in matches_minus_max:
		t = merge_tree_path(t, m, match_max)

for c in tn_uniq['c']:
	matches = search_tree(t, c, 3, [], [])
	match_len = {}
	for i in range(0, len(matches)):
		match_len[i] = len(get_from_dict(t, matches[i]))
	max_match_i = sorted(match_len.items(), key=operator.itemgetter(1), reverse=True)[0][0]
	match_max = matches[max_match_i]
	matches_minus_max = [x for x in matches if x != match_max]
	for m in matches_minus_max:
		t = merge_tree_path(t, m, match_max)

for p in tn_uniq['p']:
	matches = search_tree(t, p, 2, [], [])
	match_len = {}
	for i in range(0, len(matches)):
		match_len[i] = len(get_from_dict(t, matches[i]))
	max_match_i = sorted(match_len.items(), key=operator.itemgetter(1), reverse=True)[0][0]
	match_max = matches[max_match_i]
	matches_minus_max = [x for x in matches if x != match_max]
	for m in matches_minus_max:
		t = merge_tree_path(t, m, match_max)

for f in tn_uniq['f']:
	matches = search_tree(t, f, 5, [], [])
	match_len = {}
	for i in range(0, len(matches)):
		match_len[i] = len(get_from_dict(t, matches[i]))
	max_match_i = sorted(match_len.items(), key=operator.itemgetter(1), reverse=True)[0][0]
	match_max = matches[max_match_i]
	matches_minus_max = [x for x in matches if x != match_max]
	for m in matches_minus_max:
		t = merge_tree_path(t, m, match_max)




# output taxonomy tree cleaned
tree_out = open('canon_taxonomy_tree_clean.csv', 'w')
with tree_out:
	writer = csv.writer(tree_out)
	writer.writerow(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species list'])
	for k in t.keys():
		writer.writerow([k])
		for p in t[k].keys():
			writer.writerow(['', p])
			for c in t[k][p].keys():
				writer.writerow(['', '', c])
				for o in t[k][p][c].keys():
					writer.writerow(['', '', '', o])
					for f in t[k][p][c][o].keys():
						writer.writerow(['', '', '', '', f])
						for g in t[k][p][c][o][f].keys():
							genus_array = ['', '', '', '', '', g]
							if g in genusDict:
								genus_array.extend(genusDict[g]['sp'])
							writer.writerow(genus_array)


						






