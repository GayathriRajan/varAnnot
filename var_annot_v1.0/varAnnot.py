import requests,vcf,os,re
from operator import itemgetter
from optparse import OptionParser
from collections import OrderedDict

def get_clinvar_annot(data,chr,start,end):
    clin={'sig':'N/A'}
    print chr,start,end
    try:
        for record in data.fetch(chr,start,end):
            clin['sig']=''.join(str(e) for e in record.INFO['CLNSIG'])
        return clin
    except:
        return clin


def merge_var_annot(annot_list):
    merged = []
    for item in annot_list:
        for i, mergeditem in enumerate(merged):
            if (item[1] == mergeditem[1]) and (item[2] == mergeditem[2]):
                merged[i] += item[3:]
                break
        else:
            merged.append(item)
    return merged

def annot_var_rank(conseq):
    rank_dict = {1:['transcript_ablation','splice_acceptor_variant','splice_donor_variant'\
                'stop_gained','frameshift_variant','stop_lost','start_lost','transcript_amplification'],\
                2:['inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant'],\
                3:['splice_region_variant','incomplete_terminal_codon_variant','stop_retained_variant','synonymous_variant'],\
                4:['coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant','3_prime_UTR_variant','non_coding_transcript_exon_variant',\
                'intron_variant','NMD_transcript_variant','non_coding_transcript_variant','upstream_gene_variant','downstream_gene_variant','TFBS_ablation',\
                'TFBS_amplification','TF_binding_site_variant','regulatory_region_ablation','regulatory_region_amplification','feature_elongation',\
                'regulatory_region_variant','feature_truncation','intergenic_variant']}
    for k in rank_dict:
        if conseq in rank_dict[k]:
            yield float(str(k)+"."+str(rank_dict[k].index(conseq)))
            #yield k

def find_annot_value(d,skey,atList):
    if skey in d:
        resDict = {}
        for x in atList:
            resDict[x]=str(d.get(x).strip(":"))
            if x is 'EXON' and resDict[x]:
                resDict[x]='exon_' + resDict[x].split("/")[0]
            if x is 'INTRON' and resDict[x]:
                resDict[x]='intron_' + resDict[x].split("/")[0]
        resDict['rank'] = list(annot_var_rank(resDict['major_consequence']))[0]
        yield [resDict['rank'],resDict['major_consequence'],':'.join(resDict[x] for x in atList[1:3]),':'.join(resDict[x] for x in atList[3:])]
    for k,v in d.items():
        if isinstance(v, list):
            for i in v:
                for j in find_annot_value(i,skey,atList):
                    yield j
        elif isinstance(v, dict):
            for m in find_annot_value(v,skey,atList):
                yield m

def get_exac_csqs(url,Id):
    try:
        csqsDict = requests.get(os.path.join(url,Id)).json().get('consequence')
        csqs = [str(x) for x in csqsDict if str(x) not in '']
        #print csqs
        return csqs
    except:
        return

def get_exac_freq(url,Id):
    freqDict = {}
    features = ['rsid','exac_freq']
    for x in features:
        freqDict[x] = 'N/A'
    try:
        varFreqDict = requests.get(os.path.join(url,Id)).json()
        for x in features:
            freqDict[x] = str(varFreqDict.get(x).strip())
        #print freqDict
        return freqDict
    except:
        return freqDict

def get_exac_annot(url,Id):
    annotDict = {'type':'N/A','gene':['N/A'],'enstrans':['N/A']}
    try:
        varAnnotDict = requests.get(os.path.join(url,Id)).json()
        features = ['major_consequence','SYMBOL','HGNC_ID','HGVSc','HGVS','EXON','INTRON']
        allAnnotList = list(find_annot_value(varAnnotDict,features[0],features))
        allAnnotList = sorted(allAnnotList,key=itemgetter(0))
        allAnnotList = merge_var_annot(allAnnotList)
        annot = [elem[1:] for elem in allAnnotList if elem[0] == allAnnotList[0][0]]
        for elem in annot:
            annotDict['type'] = elem[0]
            annotDict['gene'].remove('N/A')
            annotDict['gene'].append(elem[1])
            annotDict['enstrans'].remove('N/A')
            annotDict['enstrans'].append(elem[2:])
            annotDict['gene']=[x.replace(':',' HGNC:') for x in annotDict['gene']]
            #print annotDict['gene']
        #print annotDict
        return annotDict
    except:
        return annotDict

def calc_depth_stats(rec):
    alDP = [rec['NV']] if type(rec['NV']) == int else rec['NV']
    totDP = [rec['NR']] if type(rec['NR']) == int else rec['NR']
    alRatio = ['%.2f%%' % (100.0 * a/b) for a,b in zip(alDP,totDP)]
    dpDict = {'alDP':','.join([str(x) for x in alDP]),\
           'totDP':','.join([str(x) for x in totDP]),\
           'alRatio':','.join(alRatio)}
    return dpDict

def get_var_ID(var,name):
    varDetails = var.genotype(name)
    if varDetails.is_variant == True:
        varGT = varDetails['GT'].translate(None,"|/")
        if varDetails.gt_type == 2:
            varRef=var.REF
            varAlt=var.alleles[int(varGT[0])]
        else:
            if len(varGT) == 1:
                varRef = var.REF
                varAlt = var.alleles[int(varGT)]
            else:
               (varRef,varAlt) = re.split('\||\/',varDetails.gt_bases)
    ID = '-'.join(map(str,[var.CHROM,var.POS,varRef,varAlt]))
    return ID


def main():

    parser = OptionParser(usage="usage: %prog [options]",
                          version="%prog v1.0")
    parser.add_option("-i", "--inputVCF",
                      dest="infile",
                      action="store",
                      help="input VCF file")
    
    parser.add_option("-c", "--clinvar",
                      dest="clinfile",
                      action="store",
                      help="clinvar VCF file")

    parser.add_option("-s", "--sample-name",
                      dest="samname",
                      action="store",
                      help="Sample names to analyze")
    (options, args) = parser.parse_args()

    vcfData = vcf.Reader(open(options.infile),'r')
    clinData = vcf.Reader(open(options.clinfile,'rb'))
    chr = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']
    samname = vcfData.samples
    if len(samname) == 1:
        samname = samname[0]

    indir = os.path.dirname(os.path.abspath(options.infile))
    outfile = os.path.join(indir,samname+'_annotated.txt')
    
    try:
        os.remove(outfile)
    except OSError:
        pass


    with open(outfile,'w') as rf:
        exacLink = 'http://exac.hms.harvard.edu/rest/variant/'
        rf.write("#CHROM\tPOS\tREF\tALT\trsID\tClinSig\tConseqType\tNR\tNV\tABratio\tExacAF\tGene\tTranscript"+'\n')
        for var in vcfData:
            if var.CHROM in chr:
                varDetails = var.genotype(samname)
                prtVarId = '\t'.join((map(str,[var.CHROM,var.POS,var.REF,list(var.ALT)])))
                prtVarId = prtVarId.translate(None,"[]")
                varId = get_var_ID(var,samname)
                alRatio = calc_depth_stats(varDetails)
                #print os.path.join(os.path.join(exacLink,'consequences/'),varId)
                varAnnot = get_exac_annot(os.path.join(exacLink,'consequences/'),varId)
                varFreq = get_exac_freq(os.path.join(exacLink,'variant/'),varId)
                varClin = get_clinvar_annot(clinData,var.CHROM,var.POS-1,var.POS)

                outDict = {'ID':prtVarId}
                outDict.update(alRatio)
                outDict.update(varFreq)
                outDict.update(varAnnot)
                outDict.update(varClin)
                #print outDict

                rf.write('\t'.join([outDict['ID'],outDict['rsid'],outDict['sig'],outDict['type'],outDict['alDP'],outDict['totDP'],outDict['alRatio'],outDict['exac_freq'],'\|'.join(outDict['gene']),'\|'.join([','.join(x) for x in outDict['enstrans'] if x not in ['N/A']]),'\n']))

            

main()