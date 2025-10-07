# --- wildcard constraints (no r2/clumping) ---
wildcard_constraints:
    CHR = "|".join(map(str, config["CHR"])),
    REP = "|".join(map(str, config["REP"]))

# --- final targets ---
rule all:
    input:
        # From sample_selection
        expand("{pheno_dir}/ID_filtered.txt",
               pheno_dir=config["pheno_dir"]),

        # From pheno_simulation
        expand("{pheno_dir}/simulated_pheno_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/simulated_geno_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/causal_snps{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),

        # From estimate_effects
        expand("{pheno_dir}/effects_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/training_ID_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/test_ID_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),

        # From GWAS
        expand("{results_dir}/gwas/gwas_rep{rep}_chr{chr}.pheno.glm.linear",
               results_dir=config["results_dir"],
               rep=config["REP"], chr=config["CHR"]),

        # From test_prediction
        expand("{results_dir}/test_prs_scan_rep{rep}.txt",
               results_dir=config["results_dir"], rep=config["REP"]),

        # From calculate_prs_gwas
        expand("{results_dir}/test_prs_gwas_rep{rep}.txt",
               results_dir=config["results_dir"], rep=config["REP"])


rule sample_selection:
    input:
        pheno = "/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb45105.csv",
        fam = config["fam_dir"]+"/"+"ukb22418_c1_b0_v2_s487945.fam"
    output:
        sample = config["pheno_dir"] + "/ID_filtered.txt"
    params:
        CODE=config["code_dir"]+"/"+"sample_ini_selection.R",
        N_SAMPLE = config["N_SAMPLE"],
        SAMPLE_SEED = config["SAMPLE_SEED"]
    resources:
        mem_mb = 100000,
        time = "24:00:00",
        threads = 12
    shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=12
        Rscript {params.CODE} \
                {input.fam} \
                {output.sample} \
                {params.N_SAMPLE} \
                {params.SAMPLE_SEED} \
                {input.pheno}

        module unload R/4.2.3
        """

rule filter:
    input:
        id_file = config["pheno_dir"] + "/ID_filtered.txt",
        bed = lambda wildcards: f"{config['bed_bim_dir']}/ukb22418_c{wildcards.CHR}_b0_v2_s488243.bed",
        bim = lambda wildcards: f"{config['bed_bim_dir']}/ukb22418_c{wildcards.CHR}_b0_v2_s488243.bim",
        fam = lambda wildcards: f"{config['fam_dir']}/ukb22418_c{wildcards.CHR}_b0_v2_s487945.fam"
    output:
        bed = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.bed",
        bim = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.bim",
        fam = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.fam"
    params:
        plink2 = config["plink2_path"],
        output_prefix = lambda wildcards, output: output.bed.replace('.bed', ''),
        maf = config["maf"],
        mac = config["mac"],
        geno = config["geno"],
        hwe = config["hwe"]
    threads: 8
    resources:
        mem_mb = 100000,
        time = "24:00:00",
        threads = 8
    shell:
        """
        {params.plink2} \
        --bed {input.bed} \
        --bim {input.bim} \
        --fam {input.fam} \
        --keep {input.id_file} \
        --rm-dup exclude-all \
        --snps-only just-acgt \
        --maf {params.maf} \
        --mac {params.mac} \
        --geno {params.geno} \
        --hwe {params.hwe} \
        --threads {threads} \
        --make-bed \
        --out {params.output_prefix}
        """

rule pheno_simulation:
    input:
        bim_files = expand("{pheno_dir}/ukb_chr{CHR}_filtered.bim", 
                         pheno_dir=config["pheno_dir"], CHR=config["CHR"])   
    output:
        causal_snps = config["pheno_dir"] + "/causal_snps{REP}.txt",
        pheno = config["pheno_dir"] + "/simulated_pheno_rep{REP}.txt",
        geno = config["pheno_dir"] + "/simulated_geno_rep{REP}.txt"
    params:
        REP = "{REP}",  # Wildcard for replicate number,
        seed = lambda wildcards: 1000 + int(wildcards.REP),  # Unique seed per replicate
        h2 = config["h2"],
        NUM_LOCI = config["N_LOCI"],
        plink2 = config["plink2_path"],
        NUM_ADD = config["NUM_ADDITIVE"],
        NUM_INCOM_DOM = config["NUM_INCOM_DOM"],
        NUM_COM_DOM = config["NUM_COM_DOM"],
        threads = "40",
        CODE = config["code_dir"] + "/" + "pheno_simulation_part1.R",
        pheno_dir = config["pheno_dir"],
        results_dir = config["results_dir"],
        bim_pattern = config["pheno_dir"] + "/" + "ukb_chr"
    resources:
        mem_mb = 600000,
        threads = 48,
        time = "12:00:00",
        partition = "fm-bigmem-4"
    shell:
        """
        module load R/4.2.3
        export MKL_NUM_THREADS=40
        export OMP_NUM_THREADS={params.threads}

        Rscript {params.CODE} \
                {params.bim_pattern} \
                {params.REP} \
                {params.seed} \
                {params.pheno_dir} \
                {params.results_dir} \
                {params.h2} \
                {params.NUM_LOCI} \
                {params.plink2} \
                {output.causal_snps} \
                {output.pheno} \
                {output.geno} \
                {params.NUM_ADD} \
                {params.NUM_INCOM_DOM} \
                {params.NUM_COM_DOM}
    
        module unload R/4.2.3
        """
# Test using causal SNPs only to see if the effect estimation is correct
rule estimate_effects:
    input:
        pheno = config["pheno_dir"] + "/simulated_pheno_rep{REP}.txt",
        causal_snps = config["pheno_dir"] + "/causal_snps{REP}.txt"   # ← 新增输入
    output:
        effects = config["pheno_dir"] + "/effects_rep{REP}.txt",
        train_ids = config["pheno_dir"] + "/training_ID_rep{REP}.txt",
        test_ids = config["pheno_dir"] + "/test_ID_rep{REP}.txt"
    params:
        rep = "{REP}",
        CHR = config["CHR"],
        GENO_PATTERN = config["pheno_dir"] + "/ukb_chr",
        CHUNK_SIZE = config["CHUNK_SIZE"],
        CODE = config["code_dir"] + "/" + "estimate_effects_causal.R",
        seed = lambda wildcards: 1000 + int(wildcards.REP),
        NUM_TRAIN = config["NUM_TRAIN"],
        NUM_SAMPLE = config["N_SAMPLE"],
        plink2 = config["plink2_path"]
    threads: 48
    resources:
        mem_mb = 10000,
        time = "03:00:00",
        partition = "fm-bigmem-4"
    shell:
        """
        module load R/4.2.3
        export MKL_NUM_THREADS=48
        export OMP_NUM_THREADS={threads}
      
        Rscript {params.CODE} \
                {params.rep} \
                {params.GENO_PATTERN} \
                {params.CHUNK_SIZE} \
                {output.effects} \
                {output.train_ids} \
                {output.test_ids} \
                {input.pheno} \
                {input.causal_snps} \
                {params.plink2} \
                {params.seed} \
                {params.NUM_TRAIN} \
                {params.NUM_SAMPLE}

        module unload R/4.2.3
        """

rule test_prediction:
    input:
        effects = config["pheno_dir"] + "/effects_rep{REP}.txt",
        test_ids = config["pheno_dir"] + "/test_ID_rep{REP}.txt",
        r_script = config["code_dir"] + "/test_prediction_causal.R"
    output:
        prs_scan = config["results_dir"] + "/test_prs_scan_rep{REP}.txt"
    params:
        outdir = config['results_dir'],
        CODE = config["code_dir"] + "/test_prediction_causal.R",
        GENO_PATTERN = config["pheno_dir"] + "/ukb_chr",
        plink2 = config["plink2_path"]
    threads: 20
    resources:
        mem_mb = 200000,
        time = "48:00:00"
    shell:
        """
        mkdir -p {params.outdir}
        
        module load R/4.2.3
        export MKL_NUM_THREADS={threads}
        
        Rscript {params.CODE} \
                {wildcards.REP} \
                {params.outdir} \
                {input.test_ids} \
                {input.effects} \
                {params.GENO_PATTERN} \
                {params.plink2}
            
        module unload R/4.2.3
        """


# Rule 1: create training phenotype file
rule train_pheno:
    input:
        pheno = config["pheno_dir"] + "/simulated_pheno_rep{REP}.txt",
        train_ids = config["pheno_dir"] + "/training_ID_rep{REP}.txt"
    output:
        train_pheno = config["results_dir"] + "/gwas/train_pheno_rep{REP}.txt"
    params:
        gwas_dir = config["results_dir"] + "/gwas"
    shell:
        """
        mkdir -p {params.gwas_dir}
        echo "Creating training phenotype file for rep {wildcards.REP}"
        awk 'BEGIN {{OFS="\\t"}} NR==FNR {{ids[$1]=1; next}} \
             FNR==1 {{print "FID","IID","pheno"}} \
             $2 in ids {{print $1,$2,$3}}' \
            {input.train_ids} {input.pheno} > {output.train_pheno}
        """


# Rule 2: run GWAS using prepared phenotype file
rule GWAS:
    input:
        bed = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.bed",
        bim = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.bim",
        fam = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.fam",
        train_pheno = rules.train_pheno.output.train_pheno,
        train_ids = config["pheno_dir"] + "/training_ID_rep{REP}.txt"
    output:
        gwas_out = config["results_dir"] + "/gwas/gwas_rep{REP}_chr{CHR}.pheno.glm.linear"
    params:
        out_prefix = lambda wildcards: config["results_dir"] + f"/gwas/gwas_rep{wildcards.REP}_chr{wildcards.CHR}",
        plink2 = config["plink2_path"],
        geno_prefix = lambda wildcards: config["pheno_dir"] + f"/ukb_chr{wildcards.CHR}_filtered",
        gwas_dir = config["results_dir"] + "/gwas"
    threads: 20
    resources:
        mem_mb = 100000,
        time = "24:00:00",
        partition = 'fm-bigmem-3'
    shell:
        """
        mkdir -p {params.gwas_dir}
        echo "Running GWAS for rep {wildcards.REP} chr {wildcards.CHR}"
        {params.plink2} \
            --bfile {params.geno_prefix} \
            --pheno {input.train_pheno} \
            --keep {input.train_ids} \
            --glm no-x-sex allow-no-covars 'omit-ref' \
            --no-input-missing-phenotype \
            --pheno-name pheno \
            --out {params.out_prefix} \
            --threads {threads}
        """

rule calculate_prs_gwas:
    input:
        bfile = expand("{pheno_dir}/ukb_chr{CHR}_filtered.{ext}",
                       pheno_dir=config["pheno_dir"],
                       CHR=config["CHR"],
                       ext=["bed", "bim", "fam"]),
        test_ids    = config["pheno_dir"] + "/test_ID_rep{REP}.txt",
        causal_snps = config["pheno_dir"] + "/causal_snps{REP}.txt",
        rscript     = config["code_dir"] + "/calculate_prs_gwas_causal.R",
        gwas_all    = expand(
                     config["results_dir"] + "/gwas/gwas_rep{{REP}}_chr{chr}.pheno.glm.linear",
                     chr=range(1, 23))
    output:
        prs_out = config["results_dir"] + "/test_prs_gwas_rep{REP}.txt"
    params:
        geno_pattern = config["pheno_dir"] + "/ukb_chr",
        gwas_pattern = config["results_dir"] + "/gwas/gwas_rep{REP}_chr",
        plink2       = config["plink2_path"],
        outdir       = config['results_dir'],
        rep          = "{REP}"
    threads: 40
    resources:
        mem_mb = 100000,
        time   = "24:00:00",
        partition = 'fm-bigmem-3'
    shell:
        """
        mkdir -p {params.outdir}

        module load R/4.2.3
        export MKL_NUM_THREADS={threads}

        Rscript {input.rscript} \
            "{params.geno_pattern}" \
            "{params.gwas_pattern}" \
            "{input.test_ids}" \
            "{params.plink2}" \
            "{params.outdir}" \
            "{wildcards.REP}" \
            "{output.prs_out}" \
            "{input.causal_snps}"

        module unload R/4.2.3
        """

