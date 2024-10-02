# ESM IF
## **Project Introduction**  
**ESM-IF** is a mutation assessment and inverse folding tool developed by Facebook AI Research in 2022, based on the ESM protein language model.

## Installation & Execution for HPC Deployment

### **1. Build the Environment**
```bash
conda env create -f environment.yml
conda activate inverse
conda install -c conda-forge conda-pack
conda-pack -n inverse -o inverse.tar.gz
```

### **2. Build image from docker file**
```shell
docker build -f dockerfile -t esm-if:latest .
docker build -f run_dockerfile -t run_esm-if:latest .
```

### **3. Save the Docker Image**
```shell
docker save -o run_esm-if.tar run_esm-if:latest
```

### **4. Convert Docker Image to Singularity**
```shell
singularity build run_esm-if_latest.sif docker-archive://run_esm-if.tar
```

### **5. Submit the Job to SLURM**
```shell
sbatch runhpc.sh
```