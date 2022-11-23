from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary
import pickle
import gzip
import click
import time
import logging

@click.command()
@click.argument('inputf', type=str)
@click.argument('outputf', type=str)
@click.option('--nhatm', '-n', default=75, type=int, help='number of heavy atoms')
def makelib(inputf, outputf, nhatm):
    molholder = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    patts = rdSubstructLibrary.TautomerPatternHolder()
    keys = rdSubstructLibrary.KeyFromPropHolder()
    slib = rdSubstructLibrary.SubstructLibrary(molholder,
                                                patts,
                                                keys)
    t1 = time.time()
    with gzip.GzipFile(inputf) as gz:
        suppl = Chem.ForwardSDMolSupplier(gz)
        nDone = 0
        for m in suppl:
            if m.GetNumHeavyAtoms()>nhatm:
                continue
            slib.AddMol(m)
            nDone += 1
            if nDone%5000==0 and nDone>1: print(f'did {nDone} in {time.time()-t1:.2f}s')
            
    with open(outputf, 'wb+') as outf:
        pickle.dump(slib, outf)
        print(f'That took {time.time()-t1:.2f}s in total.')

@click.command()
@click.argument('inputf', type=str)
@click.argument('ssslib', type=str)
@click.option('--nhatm', '-n', default=75, type=int, help='number of heavy atoms')
def updatelib(inputf, ssslib, nhatm=75):
    slib = pickle.load(ssslib)
    t1 = time.time()
    with gzip.GzipFile(inputf) as gz:
        suppl = Chem.ForwardSDMolSupplier(gz)
        nDone = 0
        for m in suppl:
            if m.GetNumHeavyAtoms()>nhatm:
                continue
            slib.AddMol(m)
            nDone += 1
            if nDone%5000==0 and nDone>1: print(f'did {nDone} in {time.time()-t1:.2f}s')
    with open(inputf, 'wb+') as outf:
        pickle.dump(slib, outf)
        print(f'That took {time.time()-t1:.2f}s in total.')                

@click.command()
@click.argument('qry_smarts', type=str)
@click.argument('ssslib', type=str)
@click.option('--maxresults', '-m', default=100, type=int, help='max number of results')
def runsss(qry_smarts, ssslib, maxresults):
    qry = Chem.MolFromSmarts(qry_smarts)
    slib = pickle.load(open(ssslib,'rb'))
    mids = slib.GetMatches(qry, maxResults=maxresults)
    t1 = time.time()
    with open('hits.csv', 'w') as outf:
        for mid in mids:
            outf.write(f'{Chem.MolToSmiles(slib.GetMolHolder().GetMol(mid))},{slib.GetKeyHolder().GetKey(mid)}\n')
    print(f'Done {time.time()-t1:.2f}s total')


