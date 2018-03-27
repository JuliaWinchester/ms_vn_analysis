import csv
import cPickle as pickle

class ResolvedTaxonomy():
	"""Python class for constructing MS resolved taxonomy"""
	def __init__(self, t):
		self.t = t
		self.primary_key = 1
		self.tn_id_dict = {}
		self.resolved_taxonomy = []
		self.construct_resolved_taxonomy()
	def assign_id(self, tn):
		self.tn_id_dict[tn] = self.primary_key
		self.primary_key += 1
	def tn_id(self, tn):
		if tn not in self.tn_id_dict:
			self.assign_id(tn)
		return self.tn_id_dict[tn]
	def add_tn_to_resolved_taxonomy(self, tn, rank, parent_tn):
		if parent_tn is None:
			parent_id = None
		else:
			parent_id = self.tn_id(parent_tn)
		self.resolved_taxonomy.append({
			'taxon_id': self.tn_id(tn),
			'parent_id': parent_id,
			'rank': rank,
			'name': tn
		})
	def get_nonempty_taxon(self, x):
		return next(tn for tn in reversed(x) if tn)
	def construct_resolved_taxonomy(self):
		for k in self.t.keys():
			if k:
				self.add_tn_to_resolved_taxonomy(k, 'kingdom', None)
			for p in self.t[k].keys():
				if p:
					self.add_tn_to_resolved_taxonomy(p, 'phylum', k)
				for c in self.t[k][p].keys():
					if c:
						self.add_tn_to_resolved_taxonomy(c, 'class', self.get_nonempty_taxon([k,p]))
					for o in self.t[k][p][c].keys():
						if o:
							self.add_tn_to_resolved_taxonomy(o, 'order', self.get_nonempty_taxon([k,p,c]))
						for f in self.t[k][p][c][o].keys():
							if f:
								self.add_tn_to_resolved_taxonomy(f, 'family', self.get_nonempty_taxon([k,p,c,o]))
							for g in self.t[k][p][c][o][f].keys():
								if g:
									self.add_tn_to_resolved_taxonomy(g, 'genus', self.get_nonempty_taxon([k,p,c,o,f]))
	def export_resolved_taxonomy_csv(self, filename):
		rt_out = open(filename, 'w')
		with rt_out:
			writer = csv.writer(rt_out)
			writer.writerow(['taxon_id', 'parent_id', 'rank', 'name'])
			for tn in self.resolved_taxonomy:
				writer.writerow([tn['taxon_id'], tn['parent_id'], tn['rank'], tn['name']])

t = pickle.load(open('tn_tree_manual_clean.p', 'r'))

rt = ResolvedTaxonomy(t)
rt.export_resolved_taxonomy_csv('resolved_taxonomy.csv')



								