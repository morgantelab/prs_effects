import os
from itertools import product

workdir: "/data2/morgante_lab/ukbiobank_projects/prs_effects/"
configfile: "script/replication/prs_config.yaml"

# Define wildcard constraints
wildcard_constraints:
    CHR="|".join(config["CHR"]),
    REP="|".join(config["REP"]),
    r2="|".join([f"{x:.1f}" for x in config["LD_THRESHOLDS"]]),
    p="|".join(config["P_THRESHOLDS"]),
    scenario="|".join(config["SCENARIOS"])

# Auto-create necessary directories
for scenario in config["SCENARIOS"]:
    os.makedirs(os.path.join(config["pheno_dir"], scenario), exist_ok=True)
    os.makedirs(os.path.join(config["results_dir"], scenario, "gwas"), exist_ok=True)
    for p in config["P_THRESHOLDS"]:
        os.makedirs(os.path.join(config["results_dir"], scenario, p), exist_ok=True)
        os.makedirs(os.path.join(config["results_dir"], scenario, "clumping", p), exist_ok=True)

# Entry rule that ensures the entire pipeline runs
rule all:
    input:
        expand("{pheno_dir}/{scenario}/effects_rep{rep}.txt",
               pheno_dir=config["pheno_dir"],
               scenario=config["SCENARIOS"],
               rep=config["REP"]),
        expand("{results_dir}/{scenario}/{p}/test_prs_gwas_rep{rep}.txt",
               results_dir=config["results_dir"],
               scenario=config["SCENARIOS"],
               p=config["P_THRESHOLDS"],
               rep=config["REP"]),
        expand("{results_dir}/{scenario}/gwas/gwas_rep{rep}_chr{chr}.pheno.glm.linear",
               results_dir=config["results_dir"],
               scenario=config["SCENARIOS"],
               rep=config["REP"],
               chr=config["CHR"]),
        expand("{results_dir}/{scenario}/clumping/{p}/gwas_clump_rep{rep}_chr{chr}_r2_{r2}.snplist",
               results_dir=config["results_dir"],
               scenario=config["SCENARIOS"],
               p=config["P_THRESHOLDS"],
               rep=config["REP"],
               chr=config["CHR"],
               r2=[f"{x:.1f}" for x in config["LD_THRESHOLDS"]])

rule pheno_simulation:
    input:
        expand("{pheno_dir}/ukb_chr{CHR}_filtered.bim", pheno_dir=config["pheno_dir"], CHR=config["CHR"])
    output:
        causal_snps = "{pheno_dir}/{scenario}/causal_snps{REP}.txt",
        pheno = "{pheno_dir}/{scenario}/simulated_pheno_rep{REP}.txt",
        geno = "{pheno_dir}/{scenario}/simulated_geno_rep{REP}.txt"
    params:
        seed = lambda wc: 1000 + int(wc.REP),
        NUM_LOCI = config["N_LOCI"],
        plink2 = config["plink2_path"],
        CODE = config["code_dir"] + "/pheno_simulation_part1.R",
        NUM_ADD = lambda wc: config["scenario_params"][wc.scenario]["NUM_ADDITIVE"],
        NUM_INCOM_DOM = lambda wc: config["scenario_params"][wc.scenario]["NUM_INCOM_DOM"],
        NUM_COM_DOM = lambda wc: config["scenario_params"][wc.scenario]["NUM_COM_DOM"]
    shell:
        """
        module load R/4.2.3
        Rscript {params.CODE} \
            {config[pheno_dir]}/{wildcards.scenario} \
            {wildcards.REP} \
            {params.seed} \
            {config[pheno_dir]}/{wildcards.scenario} \
            {config[results_dir]}/{wildcards.scenario} \
            {config[h2]} \
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
