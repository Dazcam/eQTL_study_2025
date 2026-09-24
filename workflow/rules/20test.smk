configfile: '../config/config.yaml' 

localrules: get_gwas_sumstats

rule all:
    input:
        expand(config["prep_gwas"]["get_gwas_sumstats"]["output"], gwas = config['gwas'])

rule get_gwas_sumstats:
    output: config["prep_gwas"]["get_gwas_sumstats"]["output"]
    params: lambda wildcards: config['gwas'][wildcards.gwas]
    message: "Download {wildcards.gwas} sumstats file"
    benchmark: "reports/benchmarks/07prep_gwas.get_gwas_sumstats_{gwas}.txt"
    log:    config["prep_gwas"]["get_gwas_sumstats"]["log"]
    run:
        import requests
        import gzip
        import os
        import io

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        with open(log[0], "w") as f_log:
            f_log.write(f"Downloading from: {params}\n")
            
            # Use a standard header to avoid being flagged as a basic script
            headers = {"User-Agent": "Mozilla/5.0"}
            
            # stream=True is vital for 200MB+ GWAS files
            with requests.get(params, stream=True, allow_redirects=True, headers=headers) as r:
                r.raise_for_status()
                
                # Check if we got JSON instead of a file
                if 'application/json' in r.headers.get('Content-Type', ''):
                    f_log.write("Error: Received JSON metadata instead of the file stream.\n")
                    raise ValueError("API returned metadata. Check if '/download' is at the end of the URL.")

                # Decompress the stream
                # We use r.raw to access the response as a file-like object
                with gzip.GzipFile(fileobj=r.raw) as gz:
                    with open(output[0], 'wb') as f_out:
                        for line in gz:
                            # Filter out '##' lines (PGC-style headers)
                            if not line.startswith(b'##'):
                                f_out.write(line)

            f_log.write("Successfully decompressed and filtered GWAS data.\n")
