dat <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/ssm_wo_cnv.txt", header=T)
dat <- dat[,-which(names(dat) == "MajorCN" | names(dat) == "MinorCN")]
dat$TotalCN <- 2
write.table(dat, "/Users/seonghwanjun/data/cell-line/bulk/OV2295/ssm_w_total_cn.txt", col.names = T, row.names = F, quote=F, sep="\t")

n_snvs <- dim(dat)[1]
dat_pwgs<-data.frame(id = paste("s", 0:(n_snvs-1), sep=""), gene = paste(dat$CHR, dat$POS, sep="_"), a = dat$d-dat$b, d = dat$d, mu_r = 0.999, mu_v = 0.499)
write.table(dat_pwgs, "/Users/seonghwanjun/data/cell-line/bulk/OV2295/PWGS/ssm.txt", quote=F, row.names=F, col.names=T, sep="\t")
id_map <- data.frame(pwgs_id=dat_pwgs$id, original_id=dat$ID)
write.table(id_map, "/Users/seonghwanjun/data/cell-line/bulk/OV2295/PWGS/id_map.txt", sep="\t", quote = F, row.names = F, col.names = T)
