'''
Converte l'AOD mergiato (ROOT) in un unico file parquet, leggendo l'albero
O2hfcandlclite da ciascuna cartella DF_* (contiene fOriginMcRec e tutte le
colonne di training) e concatenando in streaming.

run: python3 ConvertMergedAOD_to_parquet.py
'''

import uproot
import pyarrow as pa
import pyarrow.parquet as pq

# =============================================================
# CONFIGURAZIONE
# =============================================================
inputRoot   = "AO2D.part2.root"
outParquet  = "AO2D.part2.parquet"

# Solo le colonne effettivamente usate da ComputeEffScan_Lc_cent020_100pts.py
columns = [
    "fPt", "fFlagMc", "fOriginMcRec",
    "fImpactParameter0", "fImpactParameter1", "fImpactParameter2",
    "fChi2PCA", "fCpa", "fCpaXY", "fDecayLengthXY", "fDecayLength",
    "fNSigTpcTofPr0", "fNSigTpcTofPi0", "fNSigTpcTofKa1",
    "fNSigTpcTofPi2", "fNSigTpcTofPr2",
    "fM"
]

# =============================================================
# TROVA TUTTI I TREE O2hfcandlcfullp NELLE CARTELLE DF_*
# =============================================================
f = uproot.open(inputRoot)
treeKeys = [k for k in f.keys(recursive=True)
            if k.split(";")[0].endswith("/O2hfcandlclite")]

print(f"Trovati {len(treeKeys)} tree O2hfcandlclite da concatenare.")

# =============================================================
# LOOP: LEGGI OGNI DF E SCRIVI IN STREAMING SUL PARQUET
# =============================================================
writer = None
totalRows = 0

for i, key in enumerate(treeKeys):
    tree = f[key]
    arrays = tree.arrays(columns, library="ak")

    pa_table = pa.Table.from_pydict(
        {c: arrays[c].to_numpy() for c in columns}
    )

    if writer is None:
        writer = pq.ParquetWriter(outParquet, pa_table.schema)

    writer.write_table(pa_table)
    totalRows += pa_table.num_rows
    print(f"[{i+1}/{len(treeKeys)}] {key}: {pa_table.num_rows} righe scritte "
          f"(totale finora: {totalRows})")

if writer is not None:
    writer.close()

print(f"\nConversione completata: {totalRows} righe totali scritte in {outParquet}")
