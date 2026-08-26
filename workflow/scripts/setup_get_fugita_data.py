import argparse
import os
import synapseclient
import synapseutils

SYNAPSE_IDS = {
    "syn52363778": "subtypes",
    "syn52363777": "celltypes",
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--token", required=True, help="Path to file containing Synapse PAT")
    parser.add_argument("--outdir", required=True, help="Download destination directory")
    args = parser.parse_args()

    with open(args.token) as f:
        token = f.read().strip()

    syn = synapseclient.Synapse()
    syn.login(authToken=token)

    for syn_id, subdir in SYNAPSE_IDS.items():
        dest = os.path.join(args.outdir, subdir)
        os.makedirs(dest, exist_ok=True)
        print(f"Syncing {syn_id} -> {dest}")
        synapseutils.syncFromSynapse(syn, syn_id, path=dest)

if __name__ == "__main__":
    main()
