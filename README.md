# logictrn_demos
matlab code demos for developing transcription network by logictrn

LogicTRN is a quantitative framework developed to establish a mechanistic understanding on the role of combinatorial gene regulation in controlling transcriptional dynamics. This is an integrative method to decode interactions between transcription factors that form TF logics in regulating target genes. By combining logics and transcriptional kinetics into one single model framework, this tool can naturally integrate dynamic gene expression data and TFs-DNA-binding signals to identify the TF logics and to reconstruct the underlying transcriptional regulatory networks (TRNs).

### How to use
1. clone the codes from github to local disk. 
2. open the repo folder in Metlab. 
3. directely RUN the M file `logictrn_run_example.m`
4. The example files in *data* folder will be read auntomaticlly and processed by the demo code. The results of logictrn will be exported in the *output* folder. 
  + The output file will be named by `example`-`current date`-`8-digits random code`. The example output have been placed in the *output* folder (`output_example.xls`). 
  + The files in *data* folder were extracted from GSE11923 for time-series gene expression data (`gene_expression.xlsx`) and GSE39860 for DNA occupancy data (`occupancy_data.xlsx`). 
5. If you want to RUN the code with your own data, please first organised the file according to the files in *data* folder. Then CHANGE the code file `logictrn_run_example.m`, change the file path in the code.

### Need helps
If you need further helps or have any questions with this repo, please feel free to send email to the corresponding author the paper or submit to *issues*. 
