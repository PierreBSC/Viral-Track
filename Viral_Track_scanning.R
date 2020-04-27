###0)Auxiliary functions

##Function to extract information from STAR Log output

Extraction_Log_final = function(path_to_Log_file) {
  #Loading the final log file from STAR
  Log_file =read.delim(path_to_Log_file,header = F,sep="\t")
  Name_variable = as.character(Log_file$V1)
  Value_variable = as.character(Log_file$V2)
  
  
  #Extracting the informations about the mapping quality 
  Uniquely_mapped_percent = Value_variable[Name_variable=="                        Uniquely mapped reads % |"]
  Multiple_mapped_percent = Value_variable[Name_variable=="             % of reads mapped to multiple loci |"]
  Unmapped_too_short_percent = Value_variable[Name_variable=="                 % of reads unmapped: too short |"]
  Unmapped_mismatch_percent = Value_variable[Name_variable=="       % of reads unmapped: too many mismatches |"]
  Unmapped_other_percent = Value_variable[Name_variable=="                     % of reads unmapped: other |"]
  
  remove_percent = function(x) {
    l = nchar(x)
    x = substr(x,1,l-1)
    x=as.numeric(x)
    return(x)
  }  
  
  Uniquely_mapped_percent=remove_percent(Uniquely_mapped_percent)
  Multiple_mapped_percent=remove_percent(Multiple_mapped_percent)
  Unmapped_too_short_percent=remove_percent(Unmapped_too_short_percent)
  Unmapped_mismatch_percent = remove_percent(Unmapped_mismatch_percent)
  Unmapped_other_percent = remove_percent(Unmapped_other_percent)
  Total_unmapped = Unmapped_too_short_percent + Unmapped_mismatch_percent + Unmapped_other_percent
  
  Matrix_mapping = matrix(c(Uniquely_mapped_percent,Multiple_mapped_percent,Total_unmapped),nrow = 3)
  
  #####Looking at length of mappind, deletion and insertion
  
  Mean_mapped_length = as.numeric(Value_variable[Name_variable=="                          Average mapped length |"])
  Mean_deletion_length = as.numeric(Value_variable[Name_variable=="                        Deletion average length |"])
  Mean_insertion_length = as.numeric(Value_variable[Name_variable=="                       Insertion average length |"])
  #####Looking at rates of  of mismatch, insertion and deletion
  
  Mismatch_rate= remove_percent(Value_variable[Name_variable=="                      Mismatch rate per base, % |"])
  Deletion_rate = remove_percent(Value_variable[Name_variable=="                         Deletion rate per base |"])
  Insertion_rate = remove_percent(Value_variable[Name_variable=="                        Insertion rate per base |"])
  
  List_elements = list(Mapping_result= Matrix_mapping , 
                       Length_vector = c(Mean_mapped_length,Mean_insertion_length,Mean_deletion_length),
                       Rate_vector = c(Mismatch_rate,Insertion_rate,Deletion_rate) )
  
  return(List_elements)
}


##Function to compute similarity between genome sequences

alignement_score <- function(x) { # x is a vector of sequences in the Biostring format (DNAstring format)
  dist_matrix = matrix(0,nrow = length(x),ncol = length(x))
  comparison_sequence = c()
  for (i in 1:nrow(dist_matrix)){
    for (j in 1:ncol(dist_matrix)) {
      dist_matrix[i,j] = pairwiseAlignment(pattern = x[i], subject = x[j], type = "local-global", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA") , gapOpening = 1, gapExtension = 0, scoreOnly=TRUE)
    }
  }
  dist_matrix = -dist_matrix
  rownames(dist_matrix) = names(x)
  colnames(dist_matrix) = names(x)
  return(dist_matrix)
}

string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}

#Loading of the libraries

cat("Loading of the libraries.... ")

suppressMessages(library(Biostrings))
suppressMessages(library(ShortRead))
suppressMessages(library(doParallel))
suppressMessages(library(GenomicAlignments))
cat("... done ! ")


##I) First step : Parameter loading and checking

#A)Loading and pre-processing

args <- commandArgs(trailingOnly = T)

Parameter_file_path = args[1]

Parameters = read.table(Parameter_file_path,header = F,sep = "\t")
Parameters = as.character(Parameters$V1)
Parameters = strsplit(Parameters,split = "=",fixed = T)

Parameters_names = unlist(lapply(Parameters,function(x) {x[1]}))
Parameters_values = unlist(lapply(Parameters,function(x) {x[2]}))
names(Parameters_values) = Parameters_names

#B)Checking the parameters values

N_thread = as.numeric(Parameters_values["N_thread"]) #Number of threads used for the alignment by STAR
Output_directory = Parameters_values["Output_directory"] #Output directory...
Name_run = Parameters_values["Name_run"] #Name of the analytical run
Index_genome = Parameters_values["Index_genome"] #Location of the genome Index for STAR
Load_STAR_module = as.logical(Parameters_values["Load_STAR_module"]) #If STAR and samtools are stored as module : allow to load them...
Load_samtools_module = as.logical(Parameters_values["Load_samtools_module"])
Minimal_read_mapped = as.numeric(Parameters_values["Minimal_read_mapped"])
Viral_annotation_file = as.character(Parameters_values["Viral_annotation_file"])

if (is.na(N_thread)  | is.na(Output_directory) | is.na(Name_run)) {
  stop("Wrong parameter settings, check the Parameter.txt file")
}

if (!dir.exists(Output_directory)) {
  warning("Output directory does not exist ! Creating it !")
  dir.create(Output_directory)
}

if (!is.numeric(N_thread) | N_thread < 1 ) {
  stop("Wrong setting for the number of threads")
}

if (N_thread >= round(0.8*detectCores(),0) ) {
  N_thread = round(0.8*detectCores())
  stop("Not enough threads available. Check Parameter.txt")
}


if (!dir.exists(Index_genome)) {
  stop("Index genome directory does not exist")
}

if (!file.exists(Viral_annotation_file)) {
  stop("Viral database file does not exist")
}


##Registering the parallel environment
cl =makeCluster(N_thread)
registerDoParallel(cl)

##II) Second step : Target Fastq file identification

#A)Loading

Parameter_target_files = args[2]
File_to_process = read.table(Parameter_target_files,header = F,sep = "\t")
File_to_process = as.character(File_to_process$V1)
File_to_process = unique(File_to_process)

#B)Checking the values

List_target_path = c()

for (k in File_to_process) {
  if (file.exists(k)) {
    List_target_path = c(List_target_path,k)
  }
}

if (length(List_target_path)==0) {
  stop("No Fastq files identified : check the parameters !")
}

if (length(List_target_path)>0) {
  cat(paste(length(List_target_path),"Fastq files are going to be processed ! \n"))
}

##III) Third step : Mapping by itself

## If needed : adding STAR to the local environment 
if (Load_STAR_module) {
  Sys.setenv(PATH = paste("/opt/ohpc/pub/libs/star/2.7.0e/bin/Linux_x86_64_static", Sys.getenv("PATH"), sep = ":"))
}

##If neede : Adding samtools to the local environment 
if (Load_samtools_module) {
  Sys.setenv(PATH = paste("/opt/ohpc/pub/libs/samtools/1.4/bin/", Sys.getenv("PATH"), sep = ":"))
}


#For each Fastq file a sub-directory is created
#It will contains the SAM file as well as the final log file 

List_output_path = c()

for (k in List_target_path) {
  name_target = base::strsplit(k,"/",fixed = T)
  name_target = name_target[[1]]
  l = length(name_target)
  name_target = name_target[l]
  name_target = gsub('/','',name_target)
  name_target = gsub('.fastq','',name_target) #Cleaning the name to get the original Amplification batch number
  
  is_gz_file = grepl(pattern = ".gz",name_target)
  
  if (is_gz_file) {
    name_target = gsub('.gz','',name_target) #Removing if necesseray the .gz 
    
  }
  
  temp_output_dir = paste(Output_directory,"/",name_target,"/",sep = "")
  List_output_path = c(List_output_path,temp_output_dir)
  dir.create(temp_output_dir)
  cat(paste("Mapping ",name_target,".fastq file \n",sep = ""))
  
  name_prefix = paste(temp_output_dir,name_target,"_",sep = "")
  
  #We construct a complex command 
  STAR_mapping_command = paste("STAR --runThreadN",N_thread,"--genomeDir",Index_genome,"--readFilesIn",k,"--outSAMattributes NH HI AS nM NM XS ",
                               "--outFileNamePrefix",name_prefix,"--outSAMtype BAM SortedByCoordinate","--twopassMode Basic ",
                               "--outFilterMatchNmin 35 --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6")
  
  
  
  #If the file is in the format .gz then we need to add an additional paramter :
  if (is_gz_file) {
    STAR_mapping_command = paste(STAR_mapping_command, " --readFilesCommand zcat") #Allowing to read .gz files
    
  }
  
  
  system(STAR_mapping_command) ##We launch it....
  cat(paste("Mapping of",name_target,".fastq done ! \n",sep = ""))
  
}
cat("All fastq files have been mapped successfully \n")

cat("Starting the BAM file analysis \n")


##IV) Fourth step : QC and analysis of the mapping itself


Virus_database = read.delim(Viral_annotation_file,header=T,sep="\t")

for (k in List_output_path) {
  
  name_target = strsplit(k,Output_directory)
  name_target = (name_target[[1]])[2]
  name_target = gsub('/','',name_target)
  
  temp_sorted_bam = list.files(k,full.names = T,pattern = "Aligned.sortedByCoord.out.bam")[1]
  
  #To begin with : the ordered .BAM file need to indexed
  SAMtools_indexing_command = paste("samtools index",temp_sorted_bam)
  
  system(SAMtools_indexing_command)
  cat(paste("Indexing of the bam file for",name_target,"is done \n"))
  
  #Then we need to compute the number of mapped reads for each chromosome/virus
  temp_chromosome_count_path = paste(k,"/Count_chromosomes.txt",sep = "")
  
  SAMtools_chromosome_count_command = paste("samtools idxstats",temp_sorted_bam,">",temp_chromosome_count_path)
  system(SAMtools_chromosome_count_command)
  cat(paste("Computing stat file for the bam file for",name_target,"is done \n"))
  
  #We load it and clean it to remove non mapped viruses
  temp_chromosome_count = read.table(temp_chromosome_count_path,header = F,row.names = 1)
  colnames(temp_chromosome_count) = c("Chromosome_length","Mapped_reads","Unknown")
  
  ##Let's filter this table : removing host/human sequences and viruses with less than a given threshold of reads
  
  Chromosome_to_remove = c("X","Y","MT",as.character(1:23))
  Chromosome_to_remove = "1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 3, 4, 5, 6, 7, 8, 9, MT, X, Y, 
                          KI270728.1, KI270727.1, KI270442.1, KI270729.1, GL000225.1, KI270743.1, GL000008.2, GL000009.2,
                          KI270747.1, KI270722.1, GL000194.1, KI270742.1, GL000205.2, GL000195.1, KI270736.1, KI270733.1, 
                          GL000224.1, GL000219.1, KI270719.1, GL000216.2, KI270712.1, KI270706.1, KI270725.1, KI270744.1, 
                          KI270734.1, GL000213.1, GL000220.1, KI270715.1, GL000218.1, KI270749.1, KI270741.1, GL000221.1, 
                          KI270716.1, KI270731.1, KI270751.1, KI270750.1, KI270519.1, GL000214.1, KI270708.1, KI270730.1,
                          KI270438.1, KI270737.1, KI270721.1, KI270738.1, KI270748.1, KI270435.1, GL000208.1, KI270538.1,
                          KI270756.1, KI270739.1, KI270757.1, KI270709.1, KI270746.1, KI270753.1, KI270589.1, KI270726.1,
                          KI270735.1, KI270711.1, KI270745.1, KI270714.1, KI270732.1, KI270713.1, KI270754.1, KI270710.1,
                          KI270717.1, KI270724.1, KI270720.1, KI270723.1, KI270718.1, KI270317.1, KI270740.1, KI270755.1, 
                          KI270707.1, KI270579.1, KI270752.1, KI270512.1, KI270322.1, GL000226.1, KI270311.1, KI270366.1,
                          KI270511.1, KI270448.1, KI270521.1, KI270581.1, KI270582.1, KI270515.1, KI270588.1, KI270591.1, 
                          KI270522.1, KI270507.1, KI270590.1, KI270584.1, KI270320.1, KI270382.1, KI270468.1, KI270467.1,
                          KI270362.1, KI270517.1, KI270593.1, KI270528.1, KI270587.1, KI270364.1, KI270371.1, KI270333.1,
                          KI270374.1, KI270411.1, KI270414.1, KI270510.1, KI270390.1, KI270375.1, KI270420.1, KI270509.1,
                          KI270315.1, KI270302.1, KI270518.1, KI270530.1, KI270304.1, KI270418.1, KI270424.1, KI270417.1,
                          KI270508.1, KI270303.1, KI270381.1, KI270529.1, KI270425.1, KI270396.1, KI270363.1, KI270386.1,
                          KI270465.1, KI270383.1, KI270384.1, KI270330.1, KI270372.1, KI270548.1, KI270580.1, KI270387.1,
                          KI270391.1, KI270305.1, KI270373.1, KI270422.1, KI270316.1, KI270340.1, KI270338.1, KI270583.1,
                          KI270334.1, KI270429.1, KI270393.1, KI270516.1, KI270389.1, KI270466.1, KI270388.1, KI270544.1,
                          KI270310.1, KI270412.1, KI270395.1, KI270376.1, KI270337.1, KI270335.1, KI270378.1, KI270379.1,
                          KI270329.1, KI270419.1, KI270336.1, KI270312.1, KI270539.1, KI270385.1, KI270423.1, KI270392.1, KI270394.1"
  
  temp_chromosome_count = temp_chromosome_count[!rownames(temp_chromosome_count)%in%Chromosome_to_remove,] ##All viral "chromosome start with a "NC"
  temp_chromosome_count = temp_chromosome_count[temp_chromosome_count$Mapped_reads>Minimal_read_mapped,]
  
  ##We now need to check the quality of the mapping for each virus
  
  cat("Checking the mapping quality of each virus... \n")
  #We first create a sub-directory to export the sam files corresponding to each virus
  dir.create(paste(k,"Viral_BAM_files",sep = "")) 
  
  #We then create one SAM file for each virus 
  
  foreach(i=rownames(temp_chromosome_count)) %dopar% {
    if (Load_samtools_module) {
      Sys.setenv(PATH = paste("/opt/ohpc/pub/libs/samtools/1.4/bin/", Sys.getenv("PATH"), sep = ":"))
    }
    
    temp_export_bam_command = paste("samtools view -b",temp_sorted_bam,i,">",paste(k,"Viral_BAM_files/",i,".bam ",sep = ""))
    temp_export_bam_command = paste("samtools view -b ",temp_sorted_bam," \'",i,"\'"," > \'",k,"Viral_BAM_files/",i,".bam\'",sep = "")
    
    
    system(temp_export_bam_command)
    #cat(paste(rownames(temp_chromosome_count),"\n"))
  }
  cat("Export of the viral SAM file done for",name_target,"\n")
  
  #We then load them one by one and genereate a QC report 
  
  
  QC_result = foreach(i=rownames(temp_chromosome_count),.combine = rbind,.packages = c("GenomicAlignments","ShortRead")) %dopar% {
    BAM_file= readGAlignments(paste(k,"Viral_BAM_files/",i,".bam",sep = ""),param = ScanBamParam(what =scanBamWhat()))
    #Let's check the diversity of the reads
    Viral_reads = unique(BAM_file@elementMetadata$seq)
    Viral_reads_contents = alphabetFrequency(Viral_reads,as.prob =T )
    Viral_reads_contents = Viral_reads_contents[,c("A","C","G","T")]
    
    if (class(Viral_reads_contents)=="numeric") {
      Viral_reads_contents = matrix(Viral_reads_contents_mean,ncol = 4)
    }
    
    Viral_reads_contents_mean =colMeans(Viral_reads_contents)
    Read_entropy = sum(-log(Viral_reads_contents_mean)*Viral_reads_contents_mean,na.rm = T)
    
    #... the spatial distribution of the mapped reads : how much percent of the genome is mapped ?
    Covered_genome = coverage(BAM_file)[[i]]
    Covered_genome = as.numeric(Covered_genome)
    Spatial_distribution =sum(Covered_genome>0)/length(Covered_genome)
    Covered_genome = rle(sign(Covered_genome))
    Longest_contig = max(Covered_genome$lengths[Covered_genome$values>0])
    
    ##... the mean reads quality
    Reads_quality = as.character(BAM_file@elementMetadata$qual)
    Reads_quality = PhredQuality(Reads_quality)
    Reads_quality = as(Reads_quality,"IntegerList")
    Reads_quality = as.numeric(as.matrix(Reads_quality))
    Mean_read_quality = mean(Reads_quality)
    Sd_read_quality = sd(Reads_quality)
    
    ##... the number of mapped reads and unique mapped reads
    N_unique_mapped_reads = sum(BAM_file@elementMetadata$mapq==255) ##Code specific to STAR aligner.... 
    N_mapped_reads = length(BAM_file)
    Percent_uniquely_mapped = N_unique_mapped_reads/N_mapped_reads
    
    
    ##... SDUST score to compute more efficiently the quality of the reads
    Mean_dust_score = NA
    Percent_high_quality_reads = NA
    if ("ShortRead"%in%installed.packages()){
      DUST_score = dustyScore( BAM_file@elementMetadata$seq)
      Mean_dust_score = mean(DUST_score)
      Percent_high_quality_reads =  sum(DUST_score<500)/length(DUST_score)
    }
    
    ##... and lastly the pattern of the mapping 
    
    QC_temp = c(N_mapped_reads,N_unique_mapped_reads,Percent_uniquely_mapped,
                Mean_read_quality,Sd_read_quality,
                Viral_reads_contents_mean,Read_entropy,Spatial_distribution,Longest_contig,
                Mean_dust_score,Percent_high_quality_reads)
    QC_temp
  }
  
  colnames(QC_result) = c("N_reads","N_unique_reads","Percent_uniquely_mapped",
                          "Mean_read_quality","Sd_read_quality",
                          c("A","C","G","T"),"Sequence_entropy","Spatial_distribution","Longest_contig",
                          "DUST_score","Percent_high_quality_reads")
  rownames(QC_result) = rownames(temp_chromosome_count)
  QC_result = as.data.frame(QC_result)
  QC_result = QC_result[QC_result$N_unique_reads>0,]
  cat(paste("Exporting QC ",name_target,"...."))
  
  
  
  ###Now we need to extract information on the mapping by itself to get information about the QC
  
  path_to_Log_file = paste(k,name_target,"_Log.final.out",sep = "")
  Mapping_information = Extraction_Log_final(path_to_Log_file)
  Mean_mapping_length = Mapping_information$Length_vector[1]
  
  detected_virus = rownames(QC_result[QC_result$Sequence_entropy>1.2 & QC_result$Longest_contig>3*Mean_mapping_length & QC_result$Spatial_distribution>0.05,])
  
  if (length(detected_virus)==0) {
    cat("No viral sequences detected detected \n")
  }
  
  if (length(detected_virus)>0) {
    cat(paste(length(detected_virus)," viral sequences detected detected \n",sep = ""))
  }
  
  cat(paste("Exporting QC table for ",name_target,"...."))
  
  #Removing low quality virus
  
  Filtered_QC=QC_result[detected_virus,]
  
  ##Exporting the tables of the QC analysis
  write.table(file = paste(k,"QC_unfiltered.txt",sep = ""),x = QC_result,quote = F,sep = "\t")
  write.table(file = paste(k,"QC_filtered.txt",sep = ""),x = Filtered_QC,quote = F,sep = "\t")
  cat(" done !")
  
  
  
  ##But also on the splicing events identified by STAR....
  
  Splice_table_path = list.files(k,pattern = "SJ.out.tab",full.names = T)
  Splice_table = read.table(Splice_table_path,header = F,sep = "\t")
  Splice_table = Splice_table[as.character(Splice_table$V1)%in%detected_virus,]
  colnames(Splice_table) = c("Virus","Start","End","Strand","Motif","Annotated","Uniquely mapped reads","Reads crossing the junction","Unknown")
  Splice_table$Annotated[Splice_table$Annotated==1]="GT/AG"
  Splice_table$Annotated[Splice_table$Annotated==2]="CT/AC"
  Splice_table$Annotated[Splice_table$Annotated==3]="GC/AG"
  Splice_table$Annotated[Splice_table$Annotated==4]="CT/GC"
  Splice_table$Annotated[Splice_table$Annotated==5]="AT/AC"
  Splice_table$Annotated[Splice_table$Annotated==6]="GT/AT"
  Splice_table$Annotated[Splice_table$Annotated==0]="Non-canonical"
  Splice_table$Annotated = factor(Splice_table$Annotated,levels = c("GT/AG","CT/AC","GC/AG","CT/GC","AT/AC","GT/AT","Non-canonical"))
  
  ###Additional info : % of reads mapped to viral vs host
  
  Read_count_temp = read.table(temp_chromosome_count_path,header = F,row.names = 1)
  colnames(Read_count_temp) = c("Chromosome_length","Mapped_reads","Unknown")
  Read_count_temp = Read_count_temp[Read_count_temp$Mapped_reads!=0,]
  host_mapping_count = sum(Read_count_temp[grepl(pattern = "chr",rownames(Read_count_temp)),"Mapped_reads"])
  viral_mapping_count = sum(Read_count_temp[grepl(pattern = "NC",rownames(Read_count_temp)),"Mapped_reads"])
  total_mapping = viral_mapping_count + host_mapping_count
  Ratio_host_virus = matrix(data = c(host_mapping_count,viral_mapping_count)/total_mapping,ncol = 1)*100
  
  
  ##and also number of uniquely mapped reads and other reads for the filtered virus 
  
  Mapping_selected_virus = data.frame(Unique_mapping = (Filtered_QC$N_unique_reads),All_mapping = (Filtered_QC$N_reads),row.names = rownames(Filtered_QC))
  Mapping_selected_virus = Mapping_selected_virus[order(Mapping_selected_virus$Unique_mapping,decreasing = T),]
  
  
  ###Starting to plot the pdf QC
  
  cat(paste("Creating QC plot for ",name_target,"...."))
  
  
  pdf(paste(k,"QC_report.pdf",sep = ""),height = 18,width = 12)
  par(las=1,mfrow=c(4,3),mar=c(6,6,6,4))
  Color_vector = c("lightskyblue1","orange","grey80")
  #Plotting the proportion of uniquely mapped reas, unmapped etc...
  barplot(Mapping_information$Mapping_result,ylim=c(0,100),xlim=c(0,5),ylab="Percentage of reads (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5)
  legend(x = 1.5,y=50,legend = c("Unmapped","Mapped to multiple loci","Uniquely mapped"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)
  
  #Size of mapping, insertion and deletion
  barplot(Mapping_information$Length_vector,col="black",names.arg = c("Mapping length","Insertion length","Deletion length"),
          horiz = T,xlim=c(0,max(Mapping_information$Length_vector[1])*1.2),xlab="Nucleotide length",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
  
  #Rate of mismatch, deletion and insertion
  barplot(Mapping_information$Rate_vector,col="black",names.arg = c("Mismatch rate","Insertion rate","Deletion rate"),
          horiz = T,xlim=c(0,max(Mapping_information$Rate_vector[1])*1.2),xlab="Rate (%)",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
  
  #Ratio 
  Color_vector = c("darkred","grey")
  barplot(Ratio_host_virus,ylim=c(0,100),xlim=c(0,5),ylab="Mapping events (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5)
  legend(x = 1.5,y=50,legend = c("Viral mapping","Host mapping"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)
  
  #First QC for the viral hits 
  
  if (length(detected_virus) == 0) {
    Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange")) ##Viral sequences that passed QC : green
  }
  
  if (length(detected_virus) > 0) {
    Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange","green")) ##Viral sequences that passed QC : green
  }
  
  
  plot(QC_result$N_unique_reads,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,log="x",
       cex=1.5,xlab="N unique reads mapped",ylab="% Mapped genome",ylim=c(0,100),cex.lab=1.5,main="")
  abline(h=10,lwd=2,lty=2,col="grey")
  abline(v=Minimal_read_mapped,lwd=2,lty=2,col="grey")
  
  #Second QC for the viral hits 
  plot(QC_result$Sequence_entropy,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,
       cex=1.5,xlab="Sequence complexity",ylab="% Mapped genome",cex.lab=1.5,ylim=c(0,100),main="")
  abline(h=10,lwd=2,lty=2,col="grey")
  abline(v=1.2,lwd=2,lty=2,col="grey")
  
  #Third QC for the viral hits 
  
  plot(QC_result$Longest_contig,QC_result$DUST_score,pch=21,bg=Color_vector,
       cex=1.5,xlab="Longest contig (nt)",ylab="DUST score",cex.lab=1.4,main="")
  abline(v=3*Mean_mapping_length,lwd=2,lty=2,col="grey")
  
  #Number of reads for each filtered virus
  
  if (length(detected_virus) > 0) {
    
    barplot(Mapping_selected_virus$All_mapping[nrow(Mapping_selected_virus):1],
            col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
            xlim=c(0,max(Mapping_selected_virus$All_mapping)*1.2),xlab="Number of mapped reads",
            names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
    
    
    barplot(Mapping_selected_virus$Unique_mapping[nrow(Mapping_selected_virus):1],
            col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
            xlim=c(0,max(Mapping_selected_virus$Unique_mapping)*1.2),xlab="Number of uniquely mapped reads",
            names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
    
    
    #What are the identified viruses ?
    Virus_database_filtered = Virus_database[Virus_database$Name_sequence%in%detected_virus,]
    N_virus_identified = table(as.character(Virus_database_filtered$Virus_name))
    
    if (nrow(N_virus_identified)>0) {
      plot(NULL,xlim=c(0,10),ylim=c(0,length(detected_virus)),xaxt="n",yaxt="n",xlab="",ylab="")
      text(x=rep(3.5,length(detected_virus)),y=1:length(detected_virus),cex=0.8,
           Virus_database_filtered$Virus_name)
      text(x=rep(8,length(detected_virus)),y=1:length(detected_virus),cex=0.8,
           Virus_database_filtered$Name_sequence)
      
    }
    
    
  }
  
  
  dev.off()
  cat("QC plot done ! \n")
  
  
  ##Merging all viral sam files corresponding to identified viruses 
  cat("Merging Viral SAM files identified")
  
  list_BAM_files = paste(k,"/Viral_BAM_files/",sep="")
  selected_virus = list.files(list_BAM_files,full.names=F)
  selected_virus = base::strsplit(x = selected_virus,split = ".bam")
  selected_virus = unlist(lapply(selected_virus, function(x) {x[1]}))
  
  list_BAM_files = list.files(list_BAM_files,full.names=T)
  names(list_BAM_files) = selected_virus
  
  list_BAM_files = list_BAM_files[rownames(Filtered_QC)]
  list_BAM_files = paste("\'",list_BAM_files,"\'",sep = "")
  
  Merging_BAM_commad = paste("samtools merge",paste(k,"Merged_viral_mapping.bam",sep = ""),list_BAM_files)
  system(Merging_BAM_commad)
  cat("Viral detection step done !")
  
}
