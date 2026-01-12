# --------------------------------------------------------------------------------------
#
#    Gene lists and color schemes
#
# --------------------------------------------------------------------------------------

# Info  ------------------------------------------------------------------------------
# Set features
pfc_features = [
    "GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1",
    "MKI67", "C3", "ITM2A", "SST", "CALB2",
    "SCGN", "TLE3", "FEZF2", "CRYM", "LHX2"
]

wge_features = [
    "GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1",
    "MKI67", "C3", "ITM2A", "LHX6", "SIX3",
    "PROX1", "TSHZ1", "DLX1", "SCGN"
]

hip_features = [
    "NEUROD1", "GRIK4", "EOMES", "GLI3", "OLIG1",
    "MKI67", "C3", "ITM2A", "SLC17A6", "ADARB2",
    "GAD2", "TNC", "PROX1", "RELN", "LHX6"
]

tha_features = [
    "EOMES", "GLI3", "OLIG1", "MKI67", "C3",
    "ITM2A", "SLC1A6", "LHX9", "TNC", "GAD2",
    "PAX6", "SLC17A6"
]

cer_features = [
    "GAD1", "EOMES", "GLI3", "OLIG1", "MKI67",
    "C3", "ITM2A", "CA8", "ITPR1", "RBFOX3",
    "RELN"
]

general_genes = [
    'SLC17A7', 'SLC17A6', 'SLC17A8',  # VGLUT1-3
    'SLC6A1', 'SLC6A13', 'SLC6A11', 'SLC6A12',  # GABA transporters
    'SST', 'NPY', 'GAD1', 'GAD2', 'PVALB', 'CALB2', 'VIP',  # InN markers
    'C3', 'C1QB',  # MG markers
    'AQP4', 'SOX9', 'GFAP',  # Astrocytes
    'OLIG1', 'OLIG2', 'MBP',  # Oligodendrocytes
    'PDGRFA', 'PMP2',
    'EOMES', 'EBF1', 'ABCB1'
]

cholinergic_genes = [
    'CHRNA2', 'CHRNA3', 'CHRNA4', 'CHRNA5', 'CHRNA6',  # Nicotinic (brain)
    'CHRNA7', 'CHRNA9', 'CHRNA10', 'CHRNB2', 'CHRNB4',  # Nicotinic (brain)
    'CHRM1', 'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5',  # Muscarinic (brain)
    'D1BR',  # D5 Dopamine receptor (input for dopamine in Str)
    'ACHE',  # Acetylcholinesterase
    'SLC18A3',  # Vesicular acetylcholine transporter
    'ZIC4', 'LHX6', 'GBX2', 'FGF8', 'FGF17', 'DBX1'
]

claire_genes = [
    "PDGFRA",  # OPC
    "C3",  # Microglia
    "GAD1", "GAD2",  # IN
    "RELN",  # Cajal-Retzius cell
    "SLC1A3",  # Astrocytes
    "TNC", "GLI3",  # RG -> putative glioblasts
    "FN1", "FLT1",  # Endothelia
    "PAX6",  # NPC
    "MKI67",  # IP
    "RGS5",  # Pericytes
    "TJP1",  # Neuroepithelial cells
    "TLE3", "LHX2",  # upper cortical layers
    "PROX1",  # CGE
    'BCL11B',  # Sub-cortical projecting maturing Ns (Nowak 2017)
    'SATB2',  # Intercortically projecting maturing Ns (Nowak 2017)
    'CPNE8',  # Upper Layers ExNs (Nowak 2017)
    'HCRTR2',  # Subplate ExNs (Nowak 2017)
    'KCNJ6',  # Maturing ExNs (Nowak 2017)
    "RBFOX3"
]

Nowakowski_Fig4A_genes = [
    'C1QC',  # MG
    'TTR',  # Choroid Plexus
    'PECAM1',  # Endothelial Cells
    'TBX18',  # Mural Cells inc. pericytes
    'OMG',  # OPCs
    'CSPG5',  # Astrocytes
    'GHR',  # MGE-RG1
    'FZD8',  # Pallial radial glia (MGE-RG1)
    'HJURP',  # Dividing cells (RG-div1, IPC-div1, IPC-div2, MGE-div, MGE-IPC4, MGE-IPC45),
    'CRYAB',  # Truncated radial glia in G1 (tRG)
    'HEPACAM',  # Outer radial glia and astrocytes
    'LMO1',  # MGE progenitors (MGE-RG1-2, MGE-div, MGE-IPC1-3)
    'LEF1',  # Early progenitor
    'VEPH1',  # RG-like clusters among MGE progenitors (MGE-RG1-2) and CTX RG (oRG)
    'NEUROG1',  # IPC
    'NRN1',  # IPC
    'BEST3', 'PIF1', 'ASCL1', 'EOMES', 'CALB2', 'SP8', 'SST',
    'LINCO1305', 'TAC3', 'LHX6', 'NDST4', 'SLN', 'RSPO3',
    'KCNJ6', 'NEFM', 'NEFL', 'HCRTR2', 'CRYM', 'TRPM3', 'PKD1', 'UNC5D', 'NHLH1'
]

Pouliodakis_fig1c_genes = [
    'RGS5',  # Pericyte
    'CCL3', 'AIF1',  # MG
    'ITM2A', 'CLDN5', 'ESAM',  # Endo
    'SOX5',  # ExN deep layer
    'DLX1', 'DLX2', 'LHX6', 'DLX5',  # InNs
    'STMN2', 'NEUROD6', 'SATB2',  # Maturing Ns upper layer
    'PPP1R17', 'SSTR2', 'EOMES', 'PENK',  # IP
    'SATB2', 'STMN2', 'NEUROD6',  # ExN maturing
    'NEUROD6', 'POU2F2',  # Migrating ExN
    'HMGB2', 'SOX2', 'MKI67', 'PCNA',  # Cyclic Progenitors
    'PTPRZ1', 'OLIG1', 'PDGFRA',  # OPCs
    'VIM', 'PTPRZ1', 'SOX2', 'SLC1A3', 'HES1', 'HOPX'
]

# Excitatory neuron genes
exN_genes = [
    "SATB2",  # Upper-layer excitatory neurons
    "BCL11B",  # Subcortical-projecting neurons
    "TBR1",  # Deep-layer excitatory neurons
    "NEUROD6",  # Maturing excitatory neurons
    "MAP2",  # Dendritic marker
    "SYT1",  # Synaptic vesicle marker
    "GRIN2B",  # NMDA receptor subunit
    "CPNE8"  # Upper-layer maturing neurons
]

# Interneuron genes
inN_genes = [
    "GAD1",  # GABA synthesis enzyme
    "GAD2",    # GABA synthesis enzyme
    "DLX1",    # Interneuron development
    "DLX2",    # Interneuron development
    "SOX6",    # Cortical interneurons
    "CALB1",   # Calbindin
    "CALB2",   # Calretinin
    "LHX6",    # Medial ganglionic eminence interneurons
    "PROX1"    # CGE interneurons
]   

olig_genes = [
    "MBP",     # Myelin basic protein
    "MOG",     # Myelin oligodendrocyte glycoprotein
    "PLP1",    # Proteolipid protein
    "MAG",     # Myelin-associated glycoprotein
    "CNP",     # 2',3'-cyclic-nucleotide 3'-phosphodiesterase
    "OLIG1",   # Oligodendrocyte lineage
    "OLIG2"    # Oligodendrocyte lineage
]

opc_genes = [
    "PDGFRA",  # OPC marker
    "SOX10",   # Neural crest/Oligodendrocyte marker
    "NG2",     # Chondroitin sulfate proteoglycan
    "BCAN",    # Brain-specific extracellular matrix
    "TNC",     # Extracellular matrix glycoprotein
    "CSPG4"    # OPC marker
]

r_glia_genes = [
    "PAX6",    # Neural progenitor marker
    "SLC1A3",  # Astrocytic/glial glutamate transporter
    "FABP7",   # Brain fatty acid-binding protein
    "GLI3",    # Transcription factor
    "VIM",     # Intermediate filament marker
    "HES1",    # Notch signaling
    "TNC"      # Extracellular matrix glycoprotein
]

cyc_pro_genes = [
    "MKI67",   # Proliferation marker
    "TOP2A",   # DNA replication
    "PCNA",    # DNA replication and repair
    "CDK1",    # Cell cycle control
    "E2F1",    # Cell cycle regulator
    "CCNB1",   # Cyclin B1
    "AURKB"    # Aurora kinase B
]

ipc_genes = [
    "EOMES",   # Intermediate progenitors (Tbr2)
    "HOPX",    # Intermediate progenitors
    "NEUROG2", # Neurogenesis
    "NEUROD4", # Neurogenesis
    "TBR2",    # Intermediate progenitor marker
    "INSM1"    # Immature neurons/progenitors
]

mg_genes = [
    "C3",      # Complement component
    "AIF1",    # Allograft inflammatory factor 1
    "CX3CR1",  # Chemokine receptor
    "ITGAM",   # Integrin alpha M
    "P2RY12",  # Microglia-specific purinergic receptor
    "TMEM119", # Microglia marker
    "CSF1R",   # Colony-stimulating factor receptor
   "HLA-DRA"  # MHC class II
]

ExN_1 = [
    'UNC5D', 'DPY19L1', 'NRP1', 'EML6', 'DOK6', 'CLMP', 'SORBS2', 
    'CNTNAP2', 'AC119868.2', 'NRG1', 'KCNQ3'
]

ExN_3 = [
    'CUX2', 'SNTG2', 'KIF26B', 'PLXNA4', 'DLGAP2', 'TENM4', 'SLC44A5', 
    'NKAIN2', 'LIMCH1', 'FRMD4B', 'EPHA3', 'CCSER1'
]

nick_genes = [
    'PLXNA4', 'SATB2', 'CUX2',  # Upper layer ExN
    'BCL11B', 'TLE4',  # Deep layer ExN
    'GLI3',  # Radial glia
    'GAD1', 'GAD2', 'NRXN3',  # InN
    'LAMA4', 'COL4A1',  # Endothelial cells
    'PDGFRA',  # OPCs
    'C3',  # Microglia
]

final_genes = [
    "CUX2", "SATB2",            # Upper layer ExN
    "TLE4",                     # Deep layer ExN
    "GAD1", "GAD2",             # Radial glia    
    "GLI3", "PRDM16", "PAX6",   # InN
    "COL4A1", "FN1",            # Endothelial cells                  
    "PDGFRA",                   # OPCs                    
    "C3",                       # MG
    "MYT1L"                     # Extra
]

small_populations = ['VTN', 'KCNJ8', 'ABCC9', 'ART3',  # Pericytes
                     'ACTA2', 'RGS5', 'ALDH1A1',      # Smooth muscle cells
                     'SLC38A2', 'SLC4A10', 'SLC26A2', 'SLC47A1', 'FXYD5',  # Fibroblast
                     'ATP1B1', 'COL4A1', 'COL4A2', 'COL15A1', 'COLA1',
                     'COL3A1',
                     'TM4SF1', 'SLC38A5', 'CYTL1', 'BMX', 'MGP',    # Endothelial
                     'FBLN5', 'ELN', 'IGFBP4', 'CLU']

extra_from_mouse = [
    'ALDH1L1', 'SLC1A3', 'AQP4',                         # astrocyte
    'MOG', 'MAG',                                        # oligodendrocyte
    'PDGFRA', 'SUSD5', 'CSPG4',                          # OPC
    'PECAM1', 'CLDN5', 'SLCO1C1', 'OCLN',                # endothelial cell
    'GDF10', 'VIM', 'NBL1', 'A2M',                       # Bergmann glia
    'SLC17A7', 'NEUROD6', 'MAB21L1',                     # excitatory neuron
    'GAD1', 'RELN', 'CALB1',                             # inhibitory neuron
    'DES', 'MCAM', 'PDGFRB'                              # brain pericyte
]

# general
fibroblast = ['SLC38A2', 'SLC4A10', 'SLC26A2', 'SLC47A1', 'FXYD5', 
              'ATP1B1', 'COL4A1', 'COL4A2', 'COL15A1', 'COLA1',
              'COL3A1']

endothelial = ['TM4SF1', 'SLC38A5', 'CYTL1', 'BMX', 'MGP',
               'FBLN5', 'ELN', 'IGFBP4', 'CLU']

pericytes = ['VTN', 'KCNJ8', 'ABCC9', 'ART3']

bergmann = ['NPY', 'TNC', 'LINC01727', 'FST', 'MT2A', 'PIFO', 'RSPH1']
kozareva = ['PPP1R17', 'GABRA6', 'EOMES', 'LYPD6', 'PRKCD', 'SORC3', 
              'PTPRK', 'PRKCD', 'NXPH1', 'CDH22', 'KLHL1', 'ALDH1A3', 'SLC6A5', 'HTR2A', 'EDIL3',
              'DCN', 'KCNJ8', 'MRC1', 'FIT1', 'FOXJ1', 'SLC6A5', 'GRM2', 'SST', 'PTPRC']
leuko = ["PTPRC", "SKAP1", "ARHGAP15", "PRKCH", "IKZF1", "STAT4", "DOCK8", 
           "CD247", "TC2N", "IQGAP2", "FYB1", "SAMD3", "BCL11B", "CARD11", 
           "EMB", "ETS1", "HLA-E", "LCP1", "CD96", "THEMIS", "STK17B", "APBB1IP", 
           "IKZF3", "TNFAIP8", "CLEC2D", "GNG2", "CCL5", "CD53", "FLI1", 
           "ZC3HAV1"]


discreet_cols_n23 = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                       '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                       '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                       '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
                       '#cccccc', '#000000', '#ebb501']

# Colours
greens = ['#3CBB75FF', '#00FF00A5','#006400', '#B7FFB7', '#10A53DFF',
          '#95D840FF', '#9DC183',  '#708238', '#55C667FF', '#73D055FF',
          '#567D46']

ExN_blues = ['#76B5C5', '#00BDD2', '#CEE5FD', '#00B6EB', '#ABDBE3',
             '#1E81B0', '#B8D2EB', '#779CBA']

reds = ['#FAA0A0', '#FF5959', '#F75151', '#EF0029', '#D2042D']

purples = ['#B200ED',  '#DCBEFF', '#6F2DA8']


