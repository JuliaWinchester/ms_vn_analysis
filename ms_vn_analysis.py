# Script for exporting audubon core CSV file representing UF specimens.
#
# Author: Julie Winchester <julia.m.winchester@gmail.com>
# February 14, 2018

from requests import Request, Session
import credentials
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

conn = db_conn()
c = conn.cursor()

sql = """ SELECT DISTINCT genus, species FROM ms_taxonomy_names """

r = db_query(c, sql)

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






