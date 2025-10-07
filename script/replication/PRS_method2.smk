wildcard_constraints:
    CHR = "|".join(config["CHR"]),
    REP = "|".join(config["REP"]),
    r2  = "|".join([f"{x:.1f}" for x in config["LD_THRESHOLDS"]]),
    p   = "|".join(config["P_THRESHOLDS"])

import os

# 从 config 派生目录
PHENO_DIR   = config["pheno_dir"]
RESULTS_DIR = config["results_dir"]
GWAS_DIR    = config["gwas_dir"]  # 当前文件未实际使用该键，仍按 results_dir/gwas 输出

rule all:
    input:
        # Pheno directory outputs
        expand("{pheno_dir}/ID_filtered.txt",
               pheno_dir=config["pheno_dir"]),
        expand("{pheno_dir}/simulated_pheno_rep{rep}.txt", 
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/simulated_geno_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/causal_snps{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/effects_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/training_ID_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand("{pheno_dir}/test_ID_rep{rep}.txt",
               pheno_dir=config["pheno_dir"], rep=config["REP"]),
        expand(
            config["results_dir"] + "/{p}/r2_{r2}/test_prs_scan2_rep{REP}.txt",
            p=config["P_THRESHOLDS"],
            r2=[f"{x:.1f}" for x in config["LD_THRESHOLDS"]],
            REP=config["REP"]),
        # GWAS results per chromosome and replicate
        expand("{results_dir}/gwas/gwas_rep{REP}_chr{CHR}.pheno.glm.linear",
               results_dir=config["results_dir"],
               REP=config["REP"],
               CHR=config["CHR"]),
        
# －－－－－－－－－－ 删掉了原来的 rule init_dirs （不再生成目录产物） －－－－－－－－－－

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
    threads: 12
    resources:
        mem_mb = 100000,
        time = "24:00:00"
    shell:
        r"""
        mkdir -p $(dirname {output.sample})

        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS={threads}
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
        time = "24:00:00"
    shell:
        r"""
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
        REP = "{REP}",
        seed = lambda wildcards: 1000 + int(wildcards.REP),
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
        mem_mb = 200000,
        threads = 48,
        time = "12:00:00",
        partition = "fm-bigmem-4"
    shell:
        r"""
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

rule training_ID_rep_pheno:
    input:
        pheno = config["pheno_dir"] + "/simulated_pheno_rep{REP}.txt",
        train_ids = config["pheno_dir"] + "/training_ID_rep{REP}.txt"
    output:
        train_pheno = config["results_dir"] + "/gwas/train_pheno_rep{REP}.txt"
    params:
        gwas_dir = config["results_dir"] + "/gwas"
    shell:
        r"""
        mkdir -p {params.gwas_dir}
        echo "Creating training phenotype file for rep {wildcards.REP}"
        awk 'BEGIN {{OFS="\t"}} NR==FNR {{ids[$1]=1; next}} \
             FNR==1 {{print "FID","IID","pheno"}} \
             $2 in ids {{print $1,$2,$3}}' \
            {input.train_ids} {input.pheno} > {output.train_pheno}
        """

rule GWAS:
    input:
        bed = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.bed",
        bim = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.bim",
        fam = config["pheno_dir"] + "/ukb_chr{CHR}_filtered.fam",
        pheno = config["pheno_dir"] + "/simulated_pheno_rep{REP}.txt",
        # use the rule you actually defined:
        train_pheno = rules.training_ID_rep_pheno.output.train_pheno,
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
        r"""
        echo "Running GWAS for rep {wildcards.REP} chr {wildcards.CHR}"
        mkdir -p {params.gwas_dir}
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

rule union_clump_prs_scan2:
    input:
        effects   = config["pheno_dir"] + "/effects_rep{REP}.txt",
        train_ids = config["pheno_dir"] + "/training_ID_rep{REP}.txt",
        test_ids  = config["pheno_dir"] + "/test_ID_rep{REP}.txt",
        gwas = lambda wc: expand(
            f"{config['results_dir']}/gwas/gwas_rep{wc.REP}_chr{{chr}}.pheno.glm.linear",
            chr=config["CHR"]
        )
    output:
        prs_scan2 = config["results_dir"] + "/{p}/r2_{r2}/test_prs_scan2_rep{REP}.txt",
        union_tbl = config["results_dir"] + "/{p}/r2_{r2}/signif_union_rep{REP}.tsv",
        clumps    = config["results_dir"] + "/{p}/r2_{r2}/clumping_step4_beta/union_clump_rep{REP}.allchr.clumps",
        snplist   = config["results_dir"] + "/{p}/r2_{r2}/clumping_step4_beta/union_clump_rep{REP}.allchr.snplist"
    params:
        results_dir   = config["results_dir"],
        pheno_dir     = config["pheno_dir"],
        gwas_dir      = config["results_dir"] + "/gwas",
        geno_pattern  = config["pheno_dir"] + "/ukb_chr",
        plink2        = config["plink2_path"],
        clump_kb      = config.get("CLUMPING_KB", 250),
        clump_by      = "a",
        plink_mem_mb  = config.get("PLINK_MEMORY_MB", 0),
        outdir        = lambda wc: f"{config['results_dir']}/{wc.p}/r2_{wc.r2}",
        CODE          = config["code_dir"] + "/estimate_effect_v2.R"
    wildcard_constraints:
        p   = "|".join(config["P_THRESHOLDS"]),
        r2  = "|".join([f"{x:.1f}" for x in config["LD_THRESHOLDS"]]),
        REP = "|".join(config["REP"])
    threads: 20
    resources:
        mem_mb = 200000,
        time   = "48:00:00",
        partition = "fm-bigmem-3"
    shell:
        r"""
        module load R/4.2.3
        export MKL_NUM_THREADS={threads}
        Rscript {params.CODE} \
                {wildcards.REP} \
                {input.effects} \
                {input.train_ids} \
                {input.test_ids} \
                {params.gwas_dir} \
                {params.geno_pattern} \
                {params.plink2} \
                "{params.outdir}" \
                {wildcards.p} \
                {wildcards.r2} \
                {params.clump_kb} \
                {threads} \
                {params.plink_mem_mb} \
                {params.clump_by}
        module unload R/4.2.3
        """
