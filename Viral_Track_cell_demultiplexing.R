

args <- commandArgs(trailingOnly = T)


suppressMessages(library(Rsubread))
#Loading parameters
Parameter_file_path = args[1]
Parameters = read.table(Parameter_file_path,header = F,sep = "\t")
Parameters = as.character(Parameters$V1)
Parameters = strsplit(Parameters,split = "=",fixed = T)

Parameters_names = unlist(lapply(Parameters,function(x) {x[1]}))
Parameters_values = unlist(lapply(Parameters,function(x) {x[2]}))
names(Parameters_values) = Parameters_names

Output_directory = Parameters_values["Output_directory"] #Output directory
Name_run = make.names(Parameters_values["Name_run"]) #Name of the analytical run
Viral_annotation_file = as.character(Parameters_values["Viral_annotation_file"])


#Loading list files to process
Parameter_target_file = args[2]

File_to_process = read.table(Parameter_target_file,header = F,sep = "\t")
File_to_process = as.character(File_to_process$V1)
File_to_process = unique(File_to_process)

List_names = c()
for (k in File_to_process) {
  if (file.exists(k)) {
    
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
    List_names = c(List_names,name_target)
  }
}
List_target_path = paste(Output_directory,List_names,"/",sep = "/")


#Loading all viral fragments we detected in each plate

Identified_viral_fragments = c()

for (k in List_target_path) {
  name_target = names(k)
  report = read.table(paste(k,"/QC_filtered.txt",sep = ""),header = T)
  if (nrow(report)>0) {
    Identified_viral_fragments = c(Identified_viral_fragments,rownames(report))
  }
}

Identified_viral_fragments = unique(Identified_viral_fragments)
Identified_viral_fragments = Identified_viral_fragments[!is.na(Identified_viral_fragments)]


###Finding the path to the GTF : if none fonmd -> 'pseudo GTF' created

Path_GTF = paste(Output_directory,Name_run,"Merged_GTF.txt",sep = "/")

if (!file.exists(Path_GTF)) {
  cat("No GTF file found. Creating a basic GTF file ")
  
  #Loading the viral annotation file 
  Viral_annotation = read.delim(Viral_annotation_file)
  Viral_annotation = Viral_annotation[Viral_annotation$Name_sequence!=" ",]
  
  #Selecting the viral segments 
  Viral_annotation = Viral_annotation[Identified_viral_fragments,]
  Length_segments = as.numeric(Viral_annotation$Genome_length)
  
  for (k in 1:length(Identified_viral_fragments)) {
    transcript_line = c(Identified_viral_fragments[k],'Empty_Annotation','transcript',1,Length_segments[k],1000,".",".",
                        paste("gene_id ",Identified_viral_fragments[k],"_1;",""))
  }
  
}
### Counting by itself

for (path_temp in List_target_path) {
  cat("Demultiplexing reads from sample",names(List_target_path)[List_target_path==path_temp])
  
  #First aggregating all the reads from the detected viruses
  List_bam_files =c()
  
  for (segment_temp in Identified_viral_fragments) {
    List_bam_files = c(List_bam_files, paste(path_temp,"/Viral_BAM_files/",segment_temp,".bam",sep = ""))
  }
  
  command_merge = base::paste("samtools merge ",List_target_path,"/Reads_to_demultiplex.bam -f ",paste(List_bam_files,collapse = " "),sep="")
  system(command_merge)
  
  #Assigning reads to transcripts using Rsubread Featurecounts
  x = suppressMessages(featureCounts(files = paste(List_target_path,"/Reads_to_demultiplex.bam",sep=""),
                    annot.ext = Path_GTF ,isGTFAnnotationFile = T,
                    reportReads = "BAM",reportReadsPath = List_target_path,verbose = F,primaryOnly = T,allowMultiOverlap = T))
  
  #We now have to order and index the BAM file
 command_sort =paste("samtools sort ",List_target_path,"/Reads_to_demultiplex.bam.featureCounts.bam -o ",List_target_path,"/Assigned_sorted.bam",sep = "")
 system(command_sort)
 
 command_index =paste("samtools index ",List_target_path,"/Assigned_sorted.bam",sep = "")
 system(command_index)
 
 #Final command : Umi-tools command
 
 #Adding UMI_tools to the environment
 

 command_umi_tools = paste("umi_tools count --per-gene --gene-tagw=XT --assigned-status-tag=XS --per-cell -I ",
                           List_target_path,"/Assigned_sorted.bam  -S ",List_target_path, "/Expression_table.tsv  --wide-format-cell-counts",sep="")


  suppressMessages(system(command_umi_tools))
  
  
}

