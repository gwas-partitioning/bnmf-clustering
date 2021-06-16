
CURRENT <- paste(getwd(),"/",sep="")


url <- "http://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"
destfile <- paste(CURRENT,"1000GENOMES-phase_3.vcf.gz",sep="")
 
download.file(url, destfile)


library(vcfR)
vcf <- read.vcfR(destfile)

vars = vcf@fix

vars = as.data.frame(vars[,c("CHROM","POS","ID","REF","ALT")])
vars$VAR_ID = paste(vars$CHROM,vars$POS,vars$REF,vars$ALT,sep="_")
vars = vars[,c("VAR_ID","ID")]
names(vars) = c("VAR_ID","rsID")


write.table(vars,file=paste(CURRENT,"VARID_rsID_map_file.txt",sep=""), append = F, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"))
