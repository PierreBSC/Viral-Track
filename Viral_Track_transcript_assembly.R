
args <- commandArgs(trailingOnly = T)

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
names(List_target_path) = List_names

##II) Loading all viral fragments we detected in each plate

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

if (length(Identified_viral_fragments)==0) {
  stop()
  geterrmessage("No virus was detected previously. No transcript assembly is therefore possible")
}

#Creation of the Analysis directory
Path_run_directory = paste(Output_directory,Name_run,"/",sep = "/")
dir.create(Path_run_directory)

#
#Extraction 
cat("Extracting viral alignments from the different plates ... ")
for (k in Identified_viral_fragments) {
  #Creating the subdirectory for each viral sequence
  temp_path = paste(Path_run_directory,k,"/",sep = "")
  dir.create(temp_path)
  
  #Extract
  for (i in List_target_path) {
    
    name_target = gsub(Output_directory,"",i)
    name_target = gsub("/","",name_target)
    
    #Creating a large samtools command 
    samtools_extraction_command = paste("samtools view -b",paste(i,name_target,"_Aligned.sortedByCoord.out.bam",sep = "")
                                        ,k,">",paste(temp_path,name_target,".bam",sep = ""))
    system(samtools_extraction_command)
    
  }
}
cat("done !  \n")

##Merging the bam files for each viral sequence

cat("Merging the BAM files from different runs...")

for (k in Identified_viral_fragments) {
  temp_path = paste(Path_run_directory,k,"/",sep = "")
  
  samtools_merging_command = paste("samtools merge",paste(temp_path,k,"_merge.bam",sep = ""),
                                   paste(temp_path,"*",sep = ""))
  system(samtools_merging_command)
  
}

cat("...done ! \n")


cat("Performing StringTie transcriptome assembly...")


##Assembling transcriptome de novo
for (k in Identified_viral_fragments) {
  temp_path = paste(Path_run_directory,k,"/",sep = "")
  
  stringtie_command = paste("stringtie",paste(temp_path,k,"_merge.bam",sep = ""),"-f 0.01",
                            "-o",paste(temp_path,k,"_annotation.gtf",sep = ""),"-l",k)
  system(stringtie_command)
  
}

cat(" done ! \n")

cat("Merging various GTF files...")

##Merging the different annotation files into in big GTF file
Merged_GTF = c()
for (k in Identified_viral_fragments) {
  temp_path = paste(Path_run_directory,k,"/",sep = "")
  temp_path = paste(temp_path,k,"_annotation.gtf",sep = "")
  temp_GTF = try(read.delim(temp_path,skip=2,header = F))
  if(inherits(temp_GTF,"try-error")) {
    temp_GTF = NULL
    
  }
  Merged_GTF = rbind(Merged_GTF,temp_GTF)
}
write.table(Merged_GTF,paste(Path_run_directory,"/Merged_GTF.txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = F)
cat(" done ! \n")

N_transcripts_detected = sum(Merged_GTF$V3=="transcript")
cat(paste(N_transcripts_detected, "viral transcripts detected \n"))
cat("Transcript assembly step finished ! \n")


