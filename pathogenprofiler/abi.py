from Bio import SeqIO
from .utils import *
from .fasta import *
from .vcf import *
class abi:
	def __init__(self,in_obj,prefix):
		if isinstance(in_obj[0],str):
			self.filenames = in_obj
			self.records = [SeqIO.read(x,'abi') for x in self.filenames]
			self.prefix = prefix
		elif isinstance(in_obj[0],SeqIO.SeqRecord):
			self.records = in_obj
			self.prefix = prefix
		self.quals = {}
		for rec in self.records:
			self.quals[rec.name] = [ord(x) for x in rec.annotations["abif_raw"]["PCON1"]]
		self.trimmed_coords = {}
		for rec in self.records:
			trim_seq = str(SeqIO.AbiIO._abi_trim(rec).seq)
			start = str(rec.seq).index(trim_seq)
			end = start+len(trim_seq)
			self.trimmed_coords[rec.name] = list(range(start,end))
	def trim_seq(self):
		return abi([SeqIO.AbiIO._abi_trim(x) for x in self.records],self.prefix)
	def write_seq(self):
		self.fasta = "%s.fasta" % self.prefix
		with open(self.fasta,"w") as F:
			for rec in self.records:
				F.write(">%s\n%s\n" % (rec.name,rec.seq))
	def get_variants_vcf(self,refseq,gff=None):
		add_arguments_to_self(self,locals())
		self.write_seq()
		fa = fasta(self.prefix+".fasta")
		return bcf(fa.get_ref_variants(refseq,self.prefix,gff))
	def get_mixed_positions(self):
		self.channels = {'DATA9':"G", 'DATA10':"A", 'DATA11':"T", 'DATA12':"C"}
		positions = []
		for rec in self.records:
			for i,j in enumerate([rec.annotations["abif_raw"]["PLOC2"][x] for x in self.trimmed_coords[rec.name]]):
				tmp = {self.channels[x]:rec.annotations['abif_raw'][x][j] for x in self.channels}
				tot = sum(tmp.values())
				if tmp[rec.seq[i]]/tot>0.7:continue
				tmp_call = {"qual":self.quals[rec.name][i],"position":i,"alleles":tmp}
				positions.append(tmp_call)
		return positions
