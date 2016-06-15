#######01. Build the mixed library;

BML<-function(MSMSFileName="msms.txt",EvidenceFileName="evidence.txt",LibraryFileName,
              charge="no",RTShift=5,Scorediff=0,Lprob=0,Class="none"){
  ###store the names
  datamsms1<-read.table(MSMSFileName,header=T,check.names =F,stringsAsFactors=F,sep="\t")
  datamsms_name1<-names(datamsms1)
  ###
  datamsms<-read.table(MSMSFileName,header=T,stringsAsFactors=F,sep="\t")
  datamsms_name<-names(datamsms)
  
  data_evi<-read.table(EvidenceFileName,head=T,stringsAsFactors=F,sep="\t")
  #####filter msms.txt
  contamin_evi<-which(data_evi$Potential.contaminant=="+")
  contamin_pep<-unique(c(data_evi$Modified.sequence[contamin_evi]))
  ncontamin_pep<-length(contamin_pep)
  contamin_msms<-vector()
  for(ic in 1:ncontamin_pep){
    contamin_msms1<-which(datamsms$Modified.sequence==contamin_pep[ic])
    contamin_msms[ic]<-paste(contamin_msms1,collapse=";")
  }
  contamin_msms_total1<-paste(contamin_msms,collapse=";")
  contamin_msms_total<-strsplit(contamin_msms_total1,";")
  contamin_msms_total_num<-unique(c(as.numeric(contamin_msms_total[[1]])))
  ph_num<-which(datamsms$Phospho..STY.==0)
  reverse_num<-which(datamsms$Reverse=="+")
  
  num<-unique(c(ph_num,reverse_num,contamin_msms_total_num))
  newdatamsms_one<-datamsms[-num,]
  
  if(Class=="Class1"){
    score_diff_num<-which(newdatamsms_one$Score.diff>=5)
    lprob_num<-which(newdatamsms_one$Localization.prob>=0.75)
    class_index<-intersect(score_diff_num,lprob_num)
  }
  if(Class=="Class2"){
    score_diff_num<-which(newdatamsms_one$Score.diff<5)
    lprob_num<-which(newdatamsms_one$Localization.prob<0.75 & newdatamsms_one$Localization.prob>=0.25)
    class_index<-intersect(score_diff_num,lprob_num)
  }
  if(Class=="Class3"){
    score_diff_num<-which(newdatamsms_one$Score.diff<5)
    lprob_num<-which(newdatamsms_one$Localization.prob<0.25)
    class_index<-intersect(score_diff_num,lprob_num)
  }
  if(Class=="none"){
    score_diff_num<-which(newdatamsms_one$Score.diff>=Scorediff)
    lprob_num<-which(newdatamsms_one$Localization.prob>=Lprob)
    class_index<-intersect(score_diff_num,lprob_num)
  }
  
  newdatamsms<-newdatamsms_one[class_index,]
  modi_num<-which(datamsms_name=="Modified.sequence")
  newdatamsms<-newdatamsms[order(newdatamsms[,modi_num]),]
  
  #####Get the mixed library
  newdatamsms01<-newdatamsms
  if(charge=="yes"){
    newmodseq_charge<-paste(newdatamsms01$Modified.sequence,newdatamsms01$Charge,sep="-")
  }else{
    newmodseq_charge<-newdatamsms01$Modified.sequence
  }
  tabledata<-table(newmodseq_charge)
  n<-length(tabledata)
  rownum01<-vector()
  pb<-txtProgressBar(min=1,max=n,char="=>")
  for(i in 1:n){
    pms<-dimnames(tabledata)[[1]][i]
    pms_index<-which(newmodseq_charge==pms)
    npms_index<-length(pms_index)
    frt1<-vector()
    lrt<-1
    for(j in 1:npms_index){
      frt1[lrt]<-newdatamsms01$Retention.time[pms_index[j]]
      lrt<-lrt+1
    }
    rt_median1<-median(frt1)
    newfrt1<-abs(frt1-rt_median1)
    newfrt1_index<-which(newfrt1<=RTShift)
    newfrt<-frt1[newfrt1_index]
    newpms_index<-pms_index[newfrt1_index]
    nnewpms_index<-length(newpms_index)
    
    if(nnewpms_index>0){
      rownum01[i]<-newpms_index[1]
      fion<-vector()
      fintensity<-vector()
      fionmass<-vector()
      l<-1
      for(j in 1:nnewpms_index){
        fion[l]<-newdatamsms01$Matches[newpms_index[j]]
        fintensity[l]<-newdatamsms01$Intensities[newpms_index[j]]
        fionmass[l]<-newdatamsms01$Masses[newpms_index[j]]
        l<-l+1
      }
      
      newdatamsms01$Retention.time[newpms_index[1]]<-median(newfrt)
      ms_match<-paste(fion,collapse=";")
      ms_intensity<-paste(fintensity,collapse=";")
      ms_mass<-paste(fionmass,collapse=";")
      newms_match<-strsplit(ms_match,";")
      newms_intensity<-strsplit(ms_intensity,";")
      newms_mass<-strsplit(ms_mass,";")
      new_match_uni<-unique(c(newms_match[[1]]))
      newms_intensity_num<-as.numeric(newms_intensity[[1]])
      newms_mass_num<-as.numeric(newms_mass[[1]])
      nnew_match_uni<-length(new_match_uni)
      
      char_match_uni<-paste(new_match_uni,collapse=";")
      which_intensity<-vector()
      which_mass<-vector()
      for(k in 1:nnew_match_uni){
        which_newmatch<-which(newms_match[[1]]==new_match_uni[k])
        which_intensity[k]<-max(newms_intensity_num[which_newmatch])
        which_mass[k]<-newms_mass_num[which_newmatch[1]]
      }
      
      newdatamsms01$Matches[newpms_index[1]]<-char_match_uni
      newdatamsms01$Intensities[newpms_index[1]]<-paste(which_intensity,collapse=";")
      newdatamsms01$Masses[newpms_index[1]]<-paste(which_mass,collapse=";")
    }else{
      rownum01[i]<-NA
    }
    
    setTxtProgressBar(pb,i)
  }
  naomit_rownum01<-na.omit(rownum01)
  com_datamsms<-newdatamsms01[naomit_rownum01,]
  names(com_datamsms)<-datamsms_name1
  write.table(com_datamsms,file=LibraryFileName,sep="\t",row.names = F)
  close(pb)
}

########02. Calculate the isolation window according to the density of m/z;

CIW<-function(IsoWindowFileName,LibraryFileName,EvidenceFileName="evidence.txt",
              ScanrtMS1mean=0.01,ScanrtMS2mean=0.001,Peakpoints=10,
              RTSpan=c(0,120),MZSpan=c(300,1500),RTInterval=120,MZInterval=300){
  #####Get the RT length
  data_com<-read.table(LibraryFileName,head=T,stringsAsFactors=F,sep="\t")
  data_evi<-read.table(EvidenceFileName,head=T,stringsAsFactors=F,sep="\t")
  
  phnum_evi<-which(data_evi$Phospho..STY.==0)
  reverse_num_evi<-which(data_evi$Reverse=="+")
  contamin_evi<-which(data_evi$Potential.contaminant=="+")
  num_evi<-unique(c(phnum_evi,reverse_num_evi,contamin_evi))
  newdata_evi<-data_evi[-num_evi,]
  n<-length(data_com$m.z)
  rtle<-vector()
  for(ie in 1:n){
    evi_pep1<-which(newdata_evi$Modified.sequence==data_com$Modified.sequence[ie])
    if(length(evi_pep1)>0){
      rt_evi_pep1<-newdata_evi$Retention.length[evi_pep1]
      rtle[ie]<-mean(rt_evi_pep1)
      #rtle[ie]<-paste(rt_evi_pep1,collapse=";")
    }else{
      rtle[ie]<-"0"
    }
  }
  rtle_total1<-paste(rtle,collapse=";")
  rtle_total<-strsplit(rtle_total1,";")
  rtle_total_num<-as.numeric(rtle_total[[1]])
  rtle_total_num_index<-which(rtle_total_num!=0)
  rtle_total_num_nozero<-rtle_total_num[rtle_total_num_index]
  rt_quan<-mean(rtle_total_num_nozero)
  
  #####
  #####calculate circle
  npeakp<-Peakpoints
  circle_time<-rt_quan*60/npeakp
  circle_scan_num<-round((circle_time-ScanrtMS1mean*60)/(ScanrtMS2mean*60))
  #####calculate the window according to library
  isow_regulation<-1
  nrt<-RTInterval
  nmz<-MZInterval
  nnrt<-(RTSpan[2]-RTSpan[1])/nrt
  nnmz<-(MZSpan[2]-MZSpan[1])/nmz
  rt_interval<-vector()
  iso_num<-vector()
  interval_mz_num<-vector()
  interval_rt_num<-vector()
  isolationw<-vector()
  isow_center<-vector()
  for(ii in 1:nnrt){
    rt1<-0+(ii-1)*nrt
    rt2<-ii*nrt
    rt_interval[ii]<-paste(rt1,rt2,sep="_")
    rt3_num<-which(data_com$Retention.time>rt1 & data_com$Retention.time<=rt2)
    mz_total<-data_com$m.z[rt3_num]
    interval_rt_num[ii]<-sum(mz_total<=MZSpan[2])-sum(mz_total<=MZSpan[1])
    aver_total<-interval_rt_num[ii]/(MZSpan[2]-MZSpan[1])
    if(aver_total>0){
      repeat{
        seq1<-MZSpan[1]
        mz3_num<-vector()
        yyiso<-vector()
        nmzseq<-vector()
        isow_center2<-vector()
        isolationw1<-vector()
        for(jj in 1:nnmz){
          mz2<-MZSpan[1]+jj*nmz
          mz3_num[jj]<-sum(mz_total<=mz2)-sum(mz_total<=seq1)
          aver_mz3<-mz3_num[jj]/nmz
          if(aver_mz3>0){
            xx<-aver_total/aver_mz3
            yy<-round(exp(xx^(1/sqrt(xx)))+5*xx^3+isow_regulation)
            if(yy>nmz){ 
              yy<-nmz
              yyiso[jj]<-yy
              mzseq<-seq(seq1,mz2,yy)
              nmzseq[jj]<-length(mzseq)
              seq1<-mzseq[nmzseq[jj]]
            }else{
              yyiso[jj]<-yy
              mzseq<-seq(seq1,mz2,yy)
              nmzseq[jj]<-length(mzseq)
              seq1<-mzseq[nmzseq[jj]]
            }
            isow_center1<-mzseq[1:(nmzseq[jj]-1)]+yy/2
            isolationw1[jj]<-paste(mzseq,collapse=";")
            isow_center2[jj]<-paste(isow_center1,collapse=";")
          }else{
            yy<-mz2-seq1
            yyiso[jj]<-yy
            mzseq<-seq(seq1,mz2,yy)
            nmzseq[jj]<-length(mzseq)
            seq1<-mzseq[nmzseq[jj]]
            isow_center1<-mzseq[1:(nmzseq[jj]-1)]+yy/2
            isolationw1[jj]<-paste(mzseq,collapse=";")
            isow_center2[jj]<-paste(isow_center1,collapse=";")
          }
        }
        isow_regulation<-isow_regulation+1
        if((sum(nmzseq)-3)<=circle_scan_num) break
      }
      
      iso_num[ii]<-paste(yyiso,collapse=";")
      interval_mz_num[ii]<-paste(mz3_num,collapse=";")
      isolationw[ii]<-paste(isolationw1,collapse="_")
      isow_center[ii]<-paste(isow_center2,collapse="_")
    }else{
      iso_num[ii]<-"None"
      interval_mz_num[ii]<-"None"
      isolationw[ii]<-"None"
      isow_center[ii]<-"None"
    }
  }
  rt_num<-paste(interval_rt_num,collapse=";")
  isodata<-data.frame(RT=rt_interval,RT_Num=rt_num,MZ_Num=interval_mz_num,Isolation_Num=iso_num,
                      Iso_window=isolationw,Iso_window_center=isow_center)
  write.csv(isodata,file=IsoWindowFileName,row.names=F)
}

########03. Deal with the transition results from Skyline
##First, get the total peptide modified sequence that doesn't contain the "precursor"; 
##Second, extract the peptides whose n(ph) != n(STY)

ETRF<-function(TranFileName,FIonsNumber=3,qvalue=0.01){
  tran_result01<-read.csv(TranFileName,header=T,stringsAsFactors=F,sep=",")
  remove_NA01<-which(tran_result01$annotation_QValue=="#N/A")
  newtran02<-tran_result01[-remove_NA01,]
  tran_qvalue<-as.numeric(newtran02$annotation_QValue)
  Qvalue01<-which(tran_qvalue<=qvalue)
  newtran03<-newtran02[Qvalue01,]
  tran_area_index<-which(as.numeric(newtran03$Area)==0)
  newtran01<-newtran03[-tran_area_index,]
  tran_fion<-newtran01$Fragment.Ion
  ntran_fion<-length(tran_fion)
  frag_pre01<-vector()
  jfpre<-1
  for(ifpre in 1:ntran_fion){
    grep_frag_pre<-gregexpr("precursor",tran_fion[ifpre])
    if(sum(grep_frag_pre[[1]])>0){
      frag_pre01[jfpre]<-ifpre
      jfpre<-jfpre+1
    }
  }
  newtran<-newtran01[-frag_pre01,]
  
  tabledata<-table(newtran$Peptide.Modified.Sequence)
  n<-length(tabledata)
  newfname<-vector()
  newpname<-vector()
  newpms<-vector()
  newprecursor_mass<-vector()
  newcharge<-vector()
  newfion<-vector()
  newtran_mass<-vector()
  newarea<-vector()
  newrt<-vector()
  newstime<-vector()
  newetime<-vector()
  for(i in 1:n){
    pms<-dimnames(tabledata)[[1]][i]
    pms_index<-which(newtran$Peptide.Modified.Sequence==pms)
    pms_fions<-newtran$Fragment.Ion[pms_index]
    npms_fions<-length(pms_fions)
    pms_fions_new<-unique(pms_fions)
    npms_index<-length(pms_fions_new)
    
    pms_index01<-which(newtran01$Peptide.Modified.Sequence==pms)
    
    fion<-vector()
    fionmass1<-vector()
    fionmass<-vector()
    l<-1
    if(npms_index>=FIonsNumber){
      for(j in 1:npms_fions){
        fion[l]<-newtran$Fragment.Ion[pms_index[j]]
        fion_charge<-newtran$Product.Charge[pms_index[j]]
        fionmass1[l]<-newtran$Product.Mz[pms_index[j]]*fion_charge-(fion_charge-1)*1.0079
        l<-l+1
      }
      fionuni<-unique(fion)
      nfionuni<-length(fionuni)
      for(ifionuni in 1:nfionuni){
        which_fionuni<-which(fion==fionuni[ifionuni])
        fionmass[ifionuni]<-fionmass1[which_fionuni[1]]
      }
    }else{
      fionuni<-NA
      fionmass[l]<-NA
    }
    
    newfname[i]<-newtran$File.Name[pms_index[1]]
    newpname[i]<-newtran$Protein.Name[pms_index[1]]
    newpms[i]<-pms
    newprecursor_mass[i]<-newtran01$Product.Mz[pms_index01[1]]
    newcharge[i]<-newtran01$Product.Charge[pms_index01[1]]
    newfion[i]<-paste(fionuni,collapse=";")
    newtran_mass[i]<-paste(fionmass,collapse=";")
    newarea[i]<-newtran01$Total.Area.MS1[pms_index01[1]]
    newrt[i]<-newtran$Peptide.Retention.Time[pms_index[1]]
    newstime[i]<-newtran$Min.Start.Time[pms_index[1]]
    newetime[i]<-newtran$Max.End.Time[pms_index[1]]
  }
  newfion_NA<-which(newfion=="NA")
  newmsmstxt1<-data.frame(File_Name=newfname,Protein_Name=newpname,Peptide_Modified_Sequence=newpms,
                          Precursor_Mass=newprecursor_mass,Precursor_Charge=newcharge,Fragment_Ions=newfion,FIMasses=newtran_mass,
                          Area=newarea,Rentention_Time=newrt,Start_Time=newstime,End_Time=newetime)
  if(length(newfion_NA)==0){
    newmsmstxt<-newmsmstxt1
  }else{
    newmsmstxt<-newmsmstxt1[-newfion_NA,]
  }
  totfname<-paste("Outtotal_",TranFileName,sep="")
  write.csv(newmsmstxt,row.names=F,file=totfname)
  
  modseq<-as.character(newmsmstxt$Peptide_Modified_Sequence)
  modseq<-gsub("\\[\\+42\\]","a",modseq)
  modseq<-gsub("\\[\\+80\\]","p",modseq)
  modseq<-gsub("\\[\\+16\\]","o",modseq)
  modseq<-gsub("\\[\\+57\\]","c",modseq)
  modseq<-gsub("\\[\\+10\\]","r",modseq)
  modseq<-gsub("\\[\\+8\\]","l",modseq)
  modseq<-gsub("\\[.+?\\]","",modseq)
  modseqph<-gsub("a|o|c|l|r","",modseq)
  nomodseq<-gsub("a|o|c|p|l|r","",modseq)
  length_modseq<-length(newmsmstxt$Peptide_Modified_Sequence)
  
  nphsty<-vector()
  jj<-1
  for(i in 1:length_modseq){
    char_modseq_ph<-gregexpr("p",modseqph[i])
    nseqph<-sum(char_modseq_ph[[1]]>=1)
    char_sty<-gregexpr("S|T|Y",nomodseq[i])
    nchar_sty<-length(char_sty[[1]])
    no_ms_split<-strsplit(nomodseq[i],"")
    nby<-length(no_ms_split[[1]])
    char_frag_ion<-as.character(newmsmstxt$Fragment_Ions[i])
    fragion<-strsplit(char_frag_ion,";")
    
    if(nseqph==nchar_sty){
      nphsty[jj]<-i
      jj<-jj+1
    }
  }
  if(length(nphsty)>0){
    data_phsty<-newmsmstxt[-nphsty,]
  }else{
    data_phsty<-newmsmstxt
  }
  
  data_phsty_normal<-newmsmstxt[unique(nphsty),]
  normalfn<-paste("Normal_",TranFileName,sep="") #n(ph) = n(STY)
  write.csv(data_phsty_normal,row.names=F,file=normalfn)
  
  target_fn<-paste("Target_",TranFileName,sep="") #n(ph) != n(STY)
  write.csv(data_phsty,row.names=F,file=target_fn)
}

######04. Get the probable peptides list and sum the high XCorr of each ion

GPPL_Target<-function(TargetFileName,RawFileName,IsoWindowFileName,
                      min_num_fragion,CorThreshold,aboveCor,
                      ppt,peakTolThreshold){
  data_phsty_proba<-read.csv(TargetFileName,head=T,stringsAsFactors=F,sep=",")
  modseq1_first<-as.character(data_phsty_proba$Peptide_Modified_Sequence)
  modseq1_first<-gsub("\\[\\+42\\]","a",modseq1_first)
  modseq1_first<-gsub("\\[\\+80\\]","p",modseq1_first)
  modseq1_first<-gsub("\\[\\+16\\]","o",modseq1_first)
  modseq1_first<-gsub("\\[\\+57\\]","c",modseq1_first)
  modseq1_first<-gsub("\\[\\+10\\]","r",modseq1_first)
  modseq1_first<-gsub("\\[\\+8\\]","l",modseq1_first)
  modseq1_first_noph<-gsub("p","",modseq1_first)
  nomodseq1_first<-gsub("a|p|o|c|l|r","",modseq1_first)
  nmodseq1_first<-length(modseq1_first)
  File_Name<-vector()
  Protein_Name<-vector()
  Peptide_Modified_Sequence<-vector()
  Precursor_Mass<-vector()
  Precursor_Charge<-vector()
  Fragment_Ions<-vector()
  FIMasses<-vector()
  Rentention_Time<-vector()
  Start_Time<-vector()
  End_Time<-vector()
  Serial_Number_Extension<-vector()
  i_combn<-1
  pb1<- txtProgressBar(min=1,max=nmodseq1_first,char="=>")
  for(i in 1:nmodseq1_first){
    newmodseq1_first_split<-strsplit(modseq1_first_noph[i],"")
    nnewmodseq1_first<-length(newmodseq1_first_split[[1]])
    newmodseq1_first<-vector()
    for(i_first_split in 1:nnewmodseq1_first){
      if(newmodseq1_first_split[[1]][i_first_split]=="a" | newmodseq1_first_split[[1]][i_first_split]=="p" | newmodseq1_first_split[[1]][i_first_split]=="o"| newmodseq1_first_split[[1]][i_first_split]=="c"){
        newmodseq1_first[i_first_split-1]<-paste(newmodseq1_first_split[[1]][i_first_split-1],newmodseq1_first_split[[1]][i_first_split],sep="")
      }else{
        newmodseq1_first[i_first_split]<-newmodseq1_first_split[[1]][i_first_split]
      }
    }
    newmodseq2_first<-na.omit(newmodseq1_first)
    
    char_nomodseq1_first<-strsplit(nomodseq1_first[i],"")
    char_sty<-gregexpr("S|T|Y",nomodseq1_first[i])
    nchar_sty<-length(char_sty[[1]])
    char_ph<-gregexpr("p",modseq1_first[i])
    nchar_ph<-length(char_ph[[1]])
    nchoose_sty_ph<-choose(nchar_sty,nchar_ph)
    combn_sty_ph<-combinations(nchar_sty,nchar_ph)
    for(j in 1:nchoose_sty_ph){
      File_Name[i_combn]<-data_phsty_proba$File_Name[i]
      Protein_Name[i_combn]<-data_phsty_proba$Protein_Name[i]
      Precursor_Mass[i_combn]<-data_phsty_proba$Precursor_Mass[i]
      Precursor_Charge[i_combn]<-data_phsty_proba$Precursor_Charge[i]
      Fragment_Ions[i_combn]<-data_phsty_proba$Fragment_Ions[i]
      FIMasses[i_combn]<-data_phsty_proba$FIMasses[i]
      Rentention_Time[i_combn]<-data_phsty_proba$Rentention_Time[i]
      Start_Time[i_combn]<-data_phsty_proba$Start_Time[i]
      End_Time[i_combn]<-data_phsty_proba$End_Time[i]
      Serial_Number_Extension[i_combn]<-i
      
      sty_index<-char_sty[[1]][combn_sty_ph[j,]]
      new_sty<-paste(newmodseq2_first[sty_index],"p",sep="")
      newmodseq3_first<-newmodseq2_first
      newmodseq3_first[sty_index]<-new_sty
      new_char_nomodseq1_first<-paste(newmodseq3_first,collapse="")
      Peptide_Modified_Sequence[i_combn]<-new_char_nomodseq1_first
      i_combn<-i_combn+1
    }
    setTxtProgressBar(pb1,i)
  }
  new_seq1_first<-gsub("a","[+42]",Peptide_Modified_Sequence)
  new_seq1_first<-gsub("p","[+80]",new_seq1_first)
  new_seq1_first<-gsub("o","[+16]",new_seq1_first)
  new_seq1_first<-gsub("c","[+57]",new_seq1_first)
  new_seq1_first<-gsub("r","[+10]",new_seq1_first)
  new_seq1_first<-gsub("l","[+8]",new_seq1_first)
  data_phsty<-data.frame(File_Name=File_Name,
                         Protein_Name=Protein_Name,
                         Peptide_Modified_Sequence=new_seq1_first,
                         Precursor_Mass=Precursor_Mass,
                         Precursor_Charge=Precursor_Charge,
                         Fragment_Ions=Fragment_Ions,
                         FIMasses=FIMasses,
                         Rentention_Time=Rentention_Time,
                         Start_Time=Start_Time,
                         End_Time=End_Time,
                         Serial_Number_Extension=Serial_Number_Extension)
  #PSite_NoChange_FileName<-paste("Add_PSite_NoChangedFiMasses_",TargetFileName,sep="")
  #write.csv(data_phsty,row.names=F,file=PSite_NoChange_FileName)
  close(pb1)
  
  data_rawms2<-read.table(RawFileName,head=T,stringsAsFactors=F,sep="\t")
  nrawms2<-length(data_rawms2$ScanNum)
  data_isolationw<-read.csv(IsoWindowFileName,head=T,stringsAsFactors=F,sep=",")
  ndata_isolationw<-length(data_isolationw[[1]])
  niw_times<-nrawms2%/%ndata_isolationw
  newdata_isolationw<-rep(data_isolationw$Center,niw_times)
  data_isolationw_after<-data.frame(IW_Center=newdata_isolationw) #Isolation Window Center
  n_isow<-niw_times*ndata_isolationw
  newdata_rawms2<-cbind(data_rawms2[1:n_isow,],data_isolationw_after)
  
  modifiedaa<-hash(Hyd=1.0079,Car=12.0107,Nit=14.0067,Oxi=15.9994,Pho=30.9738,Sul=32.065,
                   ac=42.01056,ox=15.99491,ph=79.9663,ca=57.0215,ly=8.0142,ar=10.00827,
                   A=71.03711,R=156.10111,N=114.04293,D=115.02694,C=103.00919,
                   E=129.04259,Q=128.05858,G=57.02146,H=137.05891,I=113.08406,
                   L=113.08406,K=128.09496,M=131.04049,F=147.06841,P=97.05276,
                   S=87.03203,T=101.04768,W=186.07931,Y=163.06333,V=99.06841)
  
  n_data_phsty<-length(data_phsty$Protein_Name)
  modseq<-data_phsty$Peptide_Modified_Sequence
  modseq<-gsub("\\[\\+42\\]","a",modseq)
  modseq<-gsub("\\[\\+80\\]","p",modseq)
  modseq<-gsub("\\[\\+16\\]","o",modseq)
  modseq<-gsub("\\[\\+57\\]","c",modseq)
  modseq<-gsub("\\[\\+10\\]","r",modseq)
  modseq<-gsub("\\[\\+8\\]","l",modseq)
  modseq_ac<-gsub("p|o|c|r|l","",modseq)
  modseq_ph<-gsub("a|o|c|r|l","",modseq)
  modseq_ox<-gsub("a|p|c|r|l","",modseq)
  modseq_ca<-gsub("a|p|o|r|l","",modseq)
  modseq_ar<-gsub("a|p|o|l|c","",modseq)
  modseq_ly<-gsub("a|p|o|c|r","",modseq)
  nomodseq<-gsub("a|o|c|p|r|l","",modseq)
  pb2<-txtProgressBar(min=1,max=n_data_phsty,char="=>")
  spe_ref_XCorr_total_sum<-vector()
  spe_ref_XCorr_total_sum_filter<-vector()
  #found_frigmentions_sum_intensity<-vector()
  #found_frigmentions_sum_intensity_filter<-vector()
  found_frigmentions_sum_intensity_log<-vector()
  found_frigmentions_sum_intensity_log_filter<-vector()
  found_frigmentions_sum_intensity_total_filter<-vector()
  Number_FragIons_Total<-vector()
  number_FragIons_aboveCorthreshold<-vector()
  found_frigmentions_sum_intensity_log_filter_average<-vector()
  #WeightedCorr_total_average<-vector()
  j_proba<-1
  for(i_proba in 1:n_data_phsty){
    iw_compare1<-abs(data_phsty$Precursor_Mass[i_proba]-data_isolationw$Center[1:(ndata_isolationw/2)])
    iw_compare_index1<-which.min(iw_compare1)
    sub_ms2_iw_center_precursor1<-data_isolationw$Center[1:(ndata_isolationw/2)][iw_compare_index1]
    iw_compare2<-abs(data_phsty$Precursor_Mass[i_proba]-data_isolationw$Center[(ndata_isolationw/2+1):ndata_isolationw])
    iw_compare_index2<-which.min(iw_compare2)
    sub_ms2_iw_center_precursor2<-data_isolationw$Center[(ndata_isolationw/2+1):ndata_isolationw][iw_compare_index2]
    #for(jisow in 1:ndata_isolationw){
      #if(data_phsty$Precursor_Mass[i_proba]<data_isolationw$End[jisow] & data_phsty$Precursor_Mass[i_proba]>=data_isolationw$Start[jisow]){
        #sub_ms2_iw_center_precursor<-data_isolationw$Center[jisow]
        #break
      #}
    #}
    iw_each_index<-vector()
    iw_each_index1<-which(newdata_rawms2$IW_Center==sub_ms2_iw_center_precursor1)
    iw_each_index2<-which(newdata_rawms2$IW_Center==sub_ms2_iw_center_precursor2)
    jw<-1
    for(iw in 1:length(iw_each_index1)){
      iw_each_index[jw]<-iw_each_index1[iw]
      iw_each_index[jw+1]<-iw_each_index2[iw]
      jw<-jw+2
    }
    iw_each<-newdata_rawms2[iw_each_index,]
    datapsm01<-iw_each
    
    ####
    bmass<-vector()
    ymass<-vector()
    noms_split<-strsplit(nomodseq[i_proba],"")
    nnoms_split<-length(noms_split[[1]])
    bmass[1]<-modifiedaa[[noms_split[[1]][1]]]+modifiedaa$Hyd
    bb<-bmass[1]
    ymass[1]<-modifiedaa[[noms_split[[1]][nnoms_split]]]+modifiedaa$Hyd*3+modifiedaa$Oxi
    yy<-ymass[1]
    for(jb in 2:nnoms_split){
      bb<-bb+modifiedaa[[noms_split[[1]][jb]]]
      bmass[jb]<-bb
    }
    for(jy in 2:nnoms_split){
      yy<-yy+modifiedaa[[noms_split[[1]][nnoms_split-jy+1]]]
      ymass[jy]<-yy
    }
    
    newbmass<-bmass
    newymass<-ymass
    nby<-nnoms_split
    char_modseq_ac<-gregexpr("a",modseq_ac[i_proba])
    char_modseq_ox<-gregexpr("o",modseq_ox[i_proba])
    char_modseq_ph<-gregexpr("p",modseq_ph[i_proba])
    char_modseq_ca<-gregexpr("c",modseq_ca[i_proba])
    char_modseq_ar<-gregexpr("r",modseq_ar[i_proba])
    char_modseq_ly<-gregexpr("l",modseq_ly[i_proba])
    nseq_ac<-sum(char_modseq_ac[[1]]>=1)
    nseq_ox<-sum(char_modseq_ox[[1]]>=1)
    nseq_ph<-sum(char_modseq_ph[[1]]>=1)
    nseq_ca<-sum(char_modseq_ca[[1]]>=1)
    nseq_ar<-sum(char_modseq_ar[[1]]>=1)
    nseq_ly<-sum(char_modseq_ly[[1]]>=1)
    if(nseq_ac>=1){
      for(acj in 1:nseq_ac){
        nacj<-char_modseq_ac[[1]][acj]-acj
        newbmass[nacj:nby]<-newbmass[nacj:nby]+modifiedaa$ac
        newymass[(nby-nacj+1):nby]<-newymass[(nby-nacj+1):nby]+modifiedaa$ac
      }
    }
    if(nseq_ox>=1){
      for(oxj in 1:nseq_ox){
        noxj<-char_modseq_ox[[1]][oxj]-oxj
        newbmass[noxj:nby]<-newbmass[noxj:nby]+modifiedaa$ox
        newymass[(nby-noxj+1):nby]<-newymass[(nby-noxj+1):nby]+modifiedaa$ox
      }
    }
    if(nseq_ca>=1){
      for(caj in 1:nseq_ca){
        ncaj<-char_modseq_ca[[1]][caj]-caj
        newbmass[ncaj:nby]<-newbmass[ncaj:nby]+modifiedaa$ca
        newymass[(nby-ncaj+1):nby]<-newymass[(nby-ncaj+1):nby]+modifiedaa$ca
      }
    }
    if(nseq_ph>=1){
      for(phj in 1:nseq_ph){
        nphj<-char_modseq_ph[[1]][phj]-phj
        newbmass[nphj:nby]<-newbmass[nphj:nby]+modifiedaa$ph
        newymass[(nby-nphj+1):nby]<-newymass[(nby-nphj+1):nby]+modifiedaa$ph
      }
    }
    if(nseq_ar>=1){
      for(arj in 1:nseq_ar){
        narj<-char_modseq_ar[[1]][arj]-arj
        newbmass[narj:nby]<-newbmass[narj:nby]+modifiedaa$ar
        newymass[(nby-narj+1):nby]<-newymass[(nby-narj+1):nby]+modifiedaa$ar
      }
    }
    if(nseq_ly>=1){
      for(lyj in 1:nseq_ly){
        nlyj<-char_modseq_ly[[1]][lyj]-lyj
        newbmass[nlyj:nby]<-newbmass[nlyj:nby]+modifiedaa$ly
        newymass[(nby-nlyj+1):nby]<-newymass[(nby-nlyj+1):nby]+modifiedaa$ly
      }
    }
    newymass<-newymass-0.005
    yyrev<-rev(newymass)
    
    #### not start from the first ion
    blabel1<-paste(c("b"),1:nby,sep="")
    blabel2<-paste(blabel1," -H2O",sep="")
    blabel3<-paste(blabel1," -NH3",sep="")
    blabel4<-paste(blabel1," -98",sep="")
    
    ylabel1<-paste(c("y"),1:nby,sep="")
    ylabel2<-paste(ylabel1," -H2O",sep="")
    ylabel3<-paste(ylabel1," -NH3",sep="")
    ylabel4<-paste(ylabel1," -98",sep="")
    
    ####Obtain the decoy sequence fragment ions masses with the same b-y labels of  target peptide;
    fiby_label<-c(blabel1,blabel2,blabel3,blabel4,ylabel1,ylabel2,ylabel3,ylabel4)
    newbmass_razor<-newbmass[1:nby]
    newymass_razor<-newymass[1:nby]
    fibymass<-c(newbmass_razor,newbmass_razor-18.0152,newbmass_razor-17.0304,newbmass_razor-97.9815,
                newymass_razor,newymass_razor-18.0152,newymass_razor-17.0304,newymass_razor-97.9815)
    
    stime<-as.numeric(as.character(data_phsty$Start_Time[i_proba]))
    etime<-as.numeric(as.character(data_phsty$End_Time[i_proba]))
    newrt<-which(datapsm01$RT<=etime & datapsm01$RT>=stime)
    nnewrt<-length(newrt)
    char_sty<-gregexpr("S|T|Y",nomodseq[i_proba])
    nchar_sty<-length(char_sty[[1]])
    
    char_fi_mass<-as.character(data_phsty$FIMasses[i_proba])
    fimass_char<-strsplit(char_fi_mass,";")
    fimass<-as.numeric(fimass_char[[1]])
    nfimass<-length(fimass)
    fibymass_total<-c(fimass,fibymass)
    n_fibymass_total<-length(fibymass_total)
    
    fimass_ph1_intensity_proba<-vector()
    str_fi_intensity_proba<-vector()
    fimass_ph1_noise_proba<-vector()
    str_fi_noise_proba<-vector()
    for(sfi1_proba in 1:n_fibymass_total){
      for(nrt_proba in 1:nnewrt){
        nrt1_proba<-newrt[nrt_proba]
        char_mz_proba<-as.character(datapsm01$m.z[nrt1_proba])
        str_mz_proba<-strsplit(char_mz_proba,";")
        newstr_mz_proba<-as.numeric(str_mz_proba[[1]])
        char_intensity_proba<-as.character(datapsm01$Intensity[nrt1_proba])
        str_intensity_proba<-strsplit(char_intensity_proba,";")
        newstr_intensity_proba<-as.numeric(str_intensity_proba[[1]])
        char_noise_proba<-as.character(datapsm01$Noise[nrt1_proba])
        str_noise_proba<-strsplit(char_noise_proba,";")
        newstr_noise_proba<-as.numeric(str_noise_proba[[1]])
        
        sub_fion1_proba<-abs(newstr_mz_proba-fibymass_total[sfi1_proba])
        if(min(sub_fion1_proba)<=peakTolThreshold){
          sub_fion1_min_proba<-which.min(sub_fion1_proba)
          fimass_ph1_intensity_proba[nrt_proba]<-newstr_intensity_proba[sub_fion1_min_proba]
          fimass_ph1_noise_proba[nrt_proba]<-newstr_noise_proba[sub_fion1_min_proba]
        }else{
          fimass_ph1_intensity_proba[nrt_proba]<-0
          fimass_ph1_noise_proba[nrt_proba]<-0.0001
        }
      }
      newfimass_ph1_inten_proba<-paste(fimass_ph1_intensity_proba,collapse=";")
      str_fi_intensity_proba[sfi1_proba]<-newfimass_ph1_inten_proba
      newfimass_ph1_noise_proba<-paste(fimass_ph1_noise_proba,collapse=";")
      str_fi_noise_proba[sfi1_proba]<-newfimass_ph1_noise_proba
    }
    #get the XCorr
    spe_ref_XCorr_by_proba<-vector()
    found_frigmentions_by_intensity<-vector()
    found_frigmentions_by_intensity_log<-vector()
    ncor_by_proba<-1
    for(npby2 in (nfimass+1):n_fibymass_total){
      str_fi_intensitynew2_proba<-strsplit(str_fi_intensity_proba[npby2],";")
      str_fi_intensitynew3_proba<-as.numeric(str_fi_intensitynew2_proba[[1]])
      nnozero_str_fi_intensitynew3_1_proba<-sum(str_fi_intensitynew3_proba!=0)
      str_fi_noisenew2_proba<-strsplit(str_fi_noise_proba[npby2],";")
      str_fi_noisenew3_proba<-as.numeric(str_fi_noisenew2_proba[[1]])
      if(nnozero_str_fi_intensitynew3_1_proba<ppt){
        spe_ref_XCorr_by_proba[ncor_by_proba]<-0
        found_frigmentions_by_intensity[ncor_by_proba]<-0
        found_frigmentions_by_intensity_log[ncor_by_proba]<-0
      }else{
        cor_by1_proba<-vector()
        sum_str_fi_intensitynew1_proba<-vector()
        for(npby1_proba in 1:nfimass){
          str_fi_intensitynew_proba<-strsplit(str_fi_intensity_proba[npby1_proba],";")
          str_fi_intensitynew1_proba<-as.numeric(str_fi_intensitynew_proba[[1]])
          nnozero_str_fi_intensitynew1_proba<-sum(str_fi_intensitynew1_proba!=0)
          if(nnozero_str_fi_intensitynew1_proba<ppt){
            cor_by1_proba[npby1_proba]<-0
            sum_str_fi_intensitynew1_proba[npby1_proba]<-0
          }else{
            cor_by1_proba[npby1_proba]<-sum(str_fi_intensitynew1_proba*str_fi_intensitynew3_proba)/sqrt(sum(str_fi_intensitynew1_proba^2)*sum(str_fi_intensitynew3_proba^2))
            sum_str_fi_intensitynew1_proba[npby1_proba]<-log2(sum(str_fi_intensitynew1_proba))
          }
        }
        if(all(cor_by1_proba==0)){
          spe_ref_XCorr_by_proba[ncor_by_proba]<-0
          found_frigmentions_by_intensity[ncor_by_proba]<-0
          found_frigmentions_by_intensity_log[ncor_by_proba]<-0
        }else{
          cor_by_proba_weight<-sum(cor_by1_proba*sum_str_fi_intensitynew1_proba)/sum(sum_str_fi_intensitynew1_proba)
          cor_by_noise_proba<-str_fi_intensitynew3_proba/str_fi_noisenew3_proba
          ncor_by_noise_proba<-sum(cor_by_noise_proba>0)
          cor_by_noise_index_proba<-which(cor_by_noise_proba>0)
          cor_by_noise_sum_proba<-sum(cor_by_noise_proba[cor_by_noise_index_proba])
          spe_ref_XCorr_by_proba[ncor_by_proba]<-cor_by_proba_weight*(1-exp(-(cor_by_noise_sum_proba-ncor_by_noise_proba)))
          
          found_frigmentions_by_intensity_log[ncor_by_proba]<-log2(sum(str_fi_intensitynew3_proba))
          found_frigmentions_by_intensity[ncor_by_proba]<-sum(str_fi_intensitynew3_proba)
        }
      }
      ncor_by_proba<-ncor_by_proba+1
    }
    spe_ref_XCorr_total_sum[i_proba]<-sum(spe_ref_XCorr_by_proba)
    spe_ref_XCorr_total_sum_filter[i_proba]<-sum(spe_ref_XCorr_by_proba[which(spe_ref_XCorr_by_proba>=CorThreshold)])
    found_frigmentions_sum_intensity_log[i_proba]<-sum(found_frigmentions_by_intensity_log)
    found_frigmentions_sum_intensity_log_filter[i_proba]<-sum(found_frigmentions_by_intensity_log[which(spe_ref_XCorr_by_proba>=CorThreshold)])
    found_frigmentions_sum_intensity_total_filter[i_proba]<-sum(found_frigmentions_by_intensity[which(spe_ref_XCorr_by_proba>=CorThreshold)])
    Number_FragIons_Total[i_proba]<-sum(spe_ref_XCorr_by_proba>0)
    number_FragIons_aboveCorthreshold[i_proba]<-sum(spe_ref_XCorr_by_proba>=CorThreshold)
    found_frigmentions_sum_intensity_log_filter_average[i_proba]<-sum(found_frigmentions_by_intensity_log[which(spe_ref_XCorr_by_proba>=CorThreshold)])/sum(spe_ref_XCorr_by_proba>=CorThreshold)
    #WeightedCorr_total_average[i_proba]<-sum(spe_ref_XCorr_by_proba[which(spe_ref_XCorr_by_proba>=CorThreshold)])/sum(spe_ref_XCorr_by_proba>=CorThreshold)
    
    setTxtProgressBar(pb2,i_proba)
  }
  data_phsty_after_proba<-cbind(data_phsty,Spe_Ref_XCorr_Total_Sum=spe_ref_XCorr_total_sum,
                                Spe_Ref_XCorr_Total_Sum_Filter=spe_ref_XCorr_total_sum_filter,
                                Found_Frigmentions_Sum_Intensity_Log=found_frigmentions_sum_intensity_log,
                                Found_Frigmentions_Sum_Intensity_Log_Filter=found_frigmentions_sum_intensity_log_filter,
                                Number_FragIons_Total=Number_FragIons_Total,
                                Number_FragIons_aboveCorthreshold=number_FragIons_aboveCorthreshold,
                                Found_Frigmentions_Sum_Intensity_Log_Filter_Average=found_frigmentions_sum_intensity_log_filter_average,
                                Found_Frigmentions_Sum_Intensity_Total_Filter=found_frigmentions_sum_intensity_total_filter)
  PSite_Change_FileName<-paste("Add_PSite_ChangedFiMasses_",TargetFileName,sep="")
  write.csv(data_phsty_after_proba,row.names=F,file=PSite_Change_FileName)
  close(pb2)
  
  Serial_Number_Extension_table<-table(data_phsty_after_proba$Serial_Number_Extension)
  n_SNE_table<-length(Serial_Number_Extension_table)
  
  which_XCorr_max_paste<-vector()
  for(i_table in 1:n_SNE_table){
    XCorr_max_index<-which(data_phsty_after_proba$Serial_Number_Extension==i_table)
    XCorr_max_modpep_index<-which(data_phsty_after_proba$Peptide_Modified_Sequence[XCorr_max_index]==data_phsty_proba$Peptide_Modified_Sequence[i_table])
    modpep_WCorr_filter<-data_phsty_after_proba$Spe_Ref_XCorr_Total_Sum_Filter[XCorr_max_index][XCorr_max_modpep_index]
    modpep_intensity_average<-data_phsty_after_proba$Found_Frigmentions_Sum_Intensity_Log_Filter_Average[XCorr_max_index][XCorr_max_modpep_index]
    XCorr_max_above_index<-which(data_phsty_after_proba$Spe_Ref_XCorr_Total_Sum_Filter[XCorr_max_index]>modpep_WCorr_filter)
    if(length(XCorr_max_above_index)==0){
      which_XCorr_max_paste[i_table]<-paste(XCorr_max_index[c(XCorr_max_modpep_index)],collapse=";")
    }else{
      new_XCorr_max_above_largest_index<-XCorr_max_above_index[which.max(data_phsty_after_proba$Spe_Ref_XCorr_Total_Sum_Filter[XCorr_max_index][XCorr_max_above_index])]
      substract_XCorr<-data_phsty_after_proba$Spe_Ref_XCorr_Total_Sum_Filter[XCorr_max_index][new_XCorr_max_above_largest_index]-modpep_WCorr_filter
      if(substract_XCorr>=aboveCor){
        which_XCorr_max_paste[i_table]<-paste(XCorr_max_index[c(new_XCorr_max_above_largest_index,XCorr_max_modpep_index)],collapse=";")
      }else{
        which_XCorr_max_paste[i_table]<-paste(XCorr_max_index[c(XCorr_max_modpep_index)],collapse=";")
      }
      #new_XCorr_max_above_largest_index<-XCorr_max_above_index[which.max(data_phsty_after_proba$Spe_Ref_XCorr_Total_Sum_Filter[XCorr_max_index][XCorr_max_above_index])]
      #which_XCorr_max_paste[i_table]<-paste(XCorr_max_index[c(new_XCorr_max_above_largest_index,XCorr_max_modpep_index)],collapse=";")
    }
  }
  which_XCorr_max_paste_total<-paste(which_XCorr_max_paste,collapse=";")
  which_XCorr_max<-strsplit(which_XCorr_max_paste_total,";")
  data_phsty_after_selection<-data_phsty_after_proba[as.numeric(which_XCorr_max[[1]]),]
  data_phsty_after_selection_nozero_index<-which(data_phsty_after_selection$Number_FragIons_aboveCorthreshold>=min_num_fragion)
  new_data_phsty_after_selection<-data_phsty_after_selection[data_phsty_after_selection_nozero_index,]
  PSite_Change_FileName_selection<-paste("Add_PSite_ChangedFiMasses_selection_",TargetFileName,sep="")
  write.csv(new_data_phsty_after_selection,row.names=F,file=PSite_Change_FileName_selection)
}

GPPL_Normal<-function(NormalFileName,RawFileName,IsoWindowFileName,
                      min_num_fragion,CorThreshold,ppt,peakTolThreshold){
  data_phsty_proba<-read.csv(NormalFileName,head=T,stringsAsFactors=F,sep=",")
  modseq1_first<-as.character(data_phsty_proba$Peptide_Modified_Sequence)
  modseq1_first<-gsub("\\[\\+42\\]","a",modseq1_first)
  modseq1_first<-gsub("\\[\\+80\\]","p",modseq1_first)
  modseq1_first<-gsub("\\[\\+16\\]","o",modseq1_first)
  modseq1_first<-gsub("\\[\\+57\\]","c",modseq1_first)
  modseq1_first<-gsub("\\[\\+10\\]","r",modseq1_first)
  modseq1_first<-gsub("\\[\\+8\\]","l",modseq1_first)
  modseq1_first_noph<-gsub("p","",modseq1_first)
  nomodseq1_first<-gsub("a|p|o|c|r|l","",modseq1_first)
  nmodseq1_first<-length(modseq1_first)
  File_Name<-vector()
  Protein_Name<-vector()
  Peptide_Modified_Sequence<-vector()
  Precursor_Mass<-vector()
  Precursor_Charge<-vector()
  Fragment_Ions<-vector()
  FIMasses<-vector()
  Rentention_Time<-vector()
  Start_Time<-vector()
  End_Time<-vector()
  Serial_Number_Extension<-vector()
  i_combn<-1
  for(i in 1:nmodseq1_first){
    newmodseq1_first_split<-strsplit(modseq1_first_noph[i],"")
    nnewmodseq1_first<-length(newmodseq1_first_split[[1]])
    newmodseq1_first<-vector()
    for(i_first_split in 1:nnewmodseq1_first){
      if(newmodseq1_first_split[[1]][i_first_split]=="a" | newmodseq1_first_split[[1]][i_first_split]=="p" | newmodseq1_first_split[[1]][i_first_split]=="o"| newmodseq1_first_split[[1]][i_first_split]=="c"){
        newmodseq1_first[i_first_split-1]<-paste(newmodseq1_first_split[[1]][i_first_split-1],newmodseq1_first_split[[1]][i_first_split],sep="")
      }else{
        newmodseq1_first[i_first_split]<-newmodseq1_first_split[[1]][i_first_split]
      }
    }
    newmodseq2_first<-na.omit(newmodseq1_first)
    
    char_nomodseq1_first<-strsplit(nomodseq1_first[i],"")
    char_sty<-gregexpr("S|T|Y",nomodseq1_first[i])
    nchar_sty<-length(char_sty[[1]])
    char_ph<-gregexpr("p",modseq1_first[i])
    nchar_ph<-length(char_ph[[1]])
    nchoose_sty_ph<-choose(nchar_sty,nchar_ph)
    combn_sty_ph<-combinations(nchar_sty,nchar_ph)
    for(j in 1:nchoose_sty_ph){
      File_Name[i_combn]<-data_phsty_proba$File_Name[i]
      Protein_Name[i_combn]<-data_phsty_proba$Protein_Name[i]
      Precursor_Mass[i_combn]<-data_phsty_proba$Precursor_Mass[i]
      Precursor_Charge[i_combn]<-data_phsty_proba$Precursor_Charge[i]
      Fragment_Ions[i_combn]<-data_phsty_proba$Fragment_Ions[i]
      FIMasses[i_combn]<-data_phsty_proba$FIMasses[i]
      Rentention_Time[i_combn]<-data_phsty_proba$Rentention_Time[i]
      Start_Time[i_combn]<-data_phsty_proba$Start_Time[i]
      End_Time[i_combn]<-data_phsty_proba$End_Time[i]
      Serial_Number_Extension[i_combn]<-i
      
      sty_index<-char_sty[[1]][combn_sty_ph[j,]]
      new_sty<-paste(newmodseq2_first[sty_index],"p",sep="")
      newmodseq3_first<-newmodseq2_first
      newmodseq3_first[sty_index]<-new_sty
      new_char_nomodseq1_first<-paste(newmodseq3_first,collapse="")
      Peptide_Modified_Sequence[i_combn]<-new_char_nomodseq1_first
      i_combn<-i_combn+1
    }
  }
  new_seq1_first<-gsub("a","[+42]",Peptide_Modified_Sequence)
  new_seq1_first<-gsub("p","[+80]",new_seq1_first)
  new_seq1_first<-gsub("o","[+16]",new_seq1_first)
  new_seq1_first<-gsub("c","[+57]",new_seq1_first)
  new_seq1_first<-gsub("r","[+10]",new_seq1_first)
  new_seq1_first<-gsub("l","[+8]",new_seq1_first)
  data_phsty<-data.frame(File_Name=File_Name,
                         Protein_Name=Protein_Name,
                         Peptide_Modified_Sequence=new_seq1_first,
                         Precursor_Mass=Precursor_Mass,
                         Precursor_Charge=Precursor_Charge,
                         Fragment_Ions=Fragment_Ions,
                         FIMasses=FIMasses,
                         Rentention_Time=Rentention_Time,
                         Start_Time=Start_Time,
                         End_Time=End_Time,
                         Serial_Number_Extension=Serial_Number_Extension)
  #PSite_NoChange_FileName<-paste("Add_PSite_NoChangedFiMasses_",TargetFileName,sep="")
  #write.csv(data_phsty,row.names=F,file=PSite_NoChange_FileName)
  
  data_rawms2<-read.table(RawFileName,head=T,stringsAsFactors=F,sep="\t")
  nrawms2<-length(data_rawms2$ScanNum)
  data_isolationw<-read.csv(IsoWindowFileName,head=T,stringsAsFactors=F,sep=",")
  ndata_isolationw<-length(data_isolationw[[1]])
  niw_times<-nrawms2%/%ndata_isolationw
  newdata_isolationw<-rep(data_isolationw$Center,niw_times)
  data_isolationw_after<-data.frame(IW_Center=newdata_isolationw) #Isolation Window Center
  n_isow<-niw_times*ndata_isolationw
  newdata_rawms2<-cbind(data_rawms2[1:n_isow,],data_isolationw_after)
  
  modifiedaa<-hash(Hyd=1.0079,Car=12.0107,Nit=14.0067,Oxi=15.9994,Pho=30.9738,Sul=32.065,
                   ac=42.01056,ox=15.99491,ph=79.9663,ca=57.0215,ly=8.0142,ar=10.00827,
                   A=71.03711,R=156.10111,N=114.04293,D=115.02694,C=103.00919,
                   E=129.04259,Q=128.05858,G=57.02146,H=137.05891,I=113.08406,
                   L=113.08406,K=128.09496,M=131.04049,F=147.06841,P=97.05276,
                   S=87.03203,T=101.04768,W=186.07931,Y=163.06333,V=99.06841)
  
  n_data_phsty<-length(data_phsty$Protein_Name)
  modseq<-data_phsty$Peptide_Modified_Sequence
  modseq<-gsub("\\[\\+42\\]","a",modseq)
  modseq<-gsub("\\[\\+80\\]","p",modseq)
  modseq<-gsub("\\[\\+16\\]","o",modseq)
  modseq<-gsub("\\[\\+57\\]","c",modseq)
  modseq<-gsub("\\[\\+10\\]","r",modseq)
  modseq<-gsub("\\[\\+8\\]","l",modseq)
  modseq_ac<-gsub("p|o|c|r|l","",modseq)
  modseq_ph<-gsub("a|o|c|r|l","",modseq)
  modseq_ox<-gsub("a|p|c|r|l","",modseq)
  modseq_ca<-gsub("a|p|o|r|l","",modseq)
  modseq_ar<-gsub("a|p|o|l|c","",modseq)
  modseq_ly<-gsub("a|p|o|c|r","",modseq)
  nomodseq<-gsub("a|o|c|p|r|l","",modseq)
  #pb2<-txtProgressBar(min=1,max=n_data_phsty,char="=>")
  spe_ref_XCorr_total_sum<-vector()
  spe_ref_XCorr_total_sum_filter<-vector()
  #found_frigmentions_sum_intensity<-vector()
  #found_frigmentions_sum_intensity_filter<-vector()
  found_frigmentions_sum_intensity_log<-vector()
  found_frigmentions_sum_intensity_log_filter<-vector()
  found_frigmentions_sum_intensity_total_filter<-vector()
  Number_FragIons_Total<-vector()
  number_FragIons_aboveCorthreshold<-vector()
  found_frigmentions_sum_intensity_log_filter_average<-vector()
  #WeightedCorr_total_average<-vector()
  j_proba<-1
  for(i_proba in 1:n_data_phsty){
    iw_compare1<-abs(data_phsty$Precursor_Mass[i_proba]-data_isolationw$Center[1:(ndata_isolationw/2)])
    iw_compare_index1<-which.min(iw_compare1)
    sub_ms2_iw_center_precursor1<-data_isolationw$Center[1:(ndata_isolationw/2)][iw_compare_index1]
    iw_compare2<-abs(data_phsty$Precursor_Mass[i_proba]-data_isolationw$Center[(ndata_isolationw/2+1):ndata_isolationw])
    iw_compare_index2<-which.min(iw_compare2)
    sub_ms2_iw_center_precursor2<-data_isolationw$Center[(ndata_isolationw/2+1):ndata_isolationw][iw_compare_index2]
    #for(jisow in 1:ndata_isolationw){
    #if(data_phsty$Precursor_Mass[i_proba]<data_isolationw$End[jisow] & data_phsty$Precursor_Mass[i_proba]>=data_isolationw$Start[jisow]){
    #sub_ms2_iw_center_precursor<-data_isolationw$Center[jisow]
    #break
    #}
    #}
    iw_each_index<-vector()
    iw_each_index1<-which(newdata_rawms2$IW_Center==sub_ms2_iw_center_precursor1)
    iw_each_index2<-which(newdata_rawms2$IW_Center==sub_ms2_iw_center_precursor2)
    jw<-1
    for(iw in 1:length(iw_each_index1)){
      iw_each_index[jw]<-iw_each_index1[iw]
      iw_each_index[jw+1]<-iw_each_index2[iw]
      jw<-jw+2
    }
    iw_each<-newdata_rawms2[iw_each_index,]
    datapsm01<-iw_each
    
    ####
    bmass<-vector()
    ymass<-vector()
    noms_split<-strsplit(nomodseq[i_proba],"")
    nnoms_split<-length(noms_split[[1]])
    bmass[1]<-modifiedaa[[noms_split[[1]][1]]]+modifiedaa$Hyd
    bb<-bmass[1]
    ymass[1]<-modifiedaa[[noms_split[[1]][nnoms_split]]]+modifiedaa$Hyd*3+modifiedaa$Oxi
    yy<-ymass[1]
    for(jb in 2:nnoms_split){
      bb<-bb+modifiedaa[[noms_split[[1]][jb]]]
      bmass[jb]<-bb
    }
    for(jy in 2:nnoms_split){
      yy<-yy+modifiedaa[[noms_split[[1]][nnoms_split-jy+1]]]
      ymass[jy]<-yy
    }
    
    newbmass<-bmass
    newymass<-ymass
    nby<-nnoms_split
    char_modseq_ac<-gregexpr("a",modseq_ac[i_proba])
    char_modseq_ox<-gregexpr("o",modseq_ox[i_proba])
    char_modseq_ph<-gregexpr("p",modseq_ph[i_proba])
    char_modseq_ca<-gregexpr("c",modseq_ca[i_proba])
    char_modseq_ar<-gregexpr("r",modseq_ar[i_proba])
    char_modseq_ly<-gregexpr("l",modseq_ly[i_proba])
    nseq_ac<-sum(char_modseq_ac[[1]]>=1)
    nseq_ox<-sum(char_modseq_ox[[1]]>=1)
    nseq_ph<-sum(char_modseq_ph[[1]]>=1)
    nseq_ca<-sum(char_modseq_ca[[1]]>=1)
    nseq_ar<-sum(char_modseq_ar[[1]]>=1)
    nseq_ly<-sum(char_modseq_ly[[1]]>=1)
    if(nseq_ac>=1){
      for(acj in 1:nseq_ac){
        nacj<-char_modseq_ac[[1]][acj]-acj
        newbmass[nacj:nby]<-newbmass[nacj:nby]+modifiedaa$ac
        newymass[(nby-nacj+1):nby]<-newymass[(nby-nacj+1):nby]+modifiedaa$ac
      }
    }
    if(nseq_ox>=1){
      for(oxj in 1:nseq_ox){
        noxj<-char_modseq_ox[[1]][oxj]-oxj
        newbmass[noxj:nby]<-newbmass[noxj:nby]+modifiedaa$ox
        newymass[(nby-noxj+1):nby]<-newymass[(nby-noxj+1):nby]+modifiedaa$ox
      }
    }
    if(nseq_ca>=1){
      for(caj in 1:nseq_ca){
        ncaj<-char_modseq_ca[[1]][caj]-caj
        newbmass[ncaj:nby]<-newbmass[ncaj:nby]+modifiedaa$ca
        newymass[(nby-ncaj+1):nby]<-newymass[(nby-ncaj+1):nby]+modifiedaa$ca
      }
    }
    if(nseq_ph>=1){
      for(phj in 1:nseq_ph){
        nphj<-char_modseq_ph[[1]][phj]-phj
        newbmass[nphj:nby]<-newbmass[nphj:nby]+modifiedaa$ph
        newymass[(nby-nphj+1):nby]<-newymass[(nby-nphj+1):nby]+modifiedaa$ph
      }
    }
    if(nseq_ar>=1){
      for(arj in 1:nseq_ar){
        narj<-char_modseq_ar[[1]][arj]-arj
        newbmass[narj:nby]<-newbmass[narj:nby]+modifiedaa$ar
        newymass[(nby-narj+1):nby]<-newymass[(nby-narj+1):nby]+modifiedaa$ar
      }
    }
    if(nseq_ly>=1){
      for(lyj in 1:nseq_ly){
        nlyj<-char_modseq_ly[[1]][lyj]-lyj
        newbmass[nlyj:nby]<-newbmass[nlyj:nby]+modifiedaa$ly
        newymass[(nby-nlyj+1):nby]<-newymass[(nby-nlyj+1):nby]+modifiedaa$ly
      }
    }
    newymass<-newymass-0.005
    yyrev<-rev(newymass)
    
    #### not start from the first ion
    blabel1<-paste(c("b"),1:nby,sep="")
    blabel2<-paste(blabel1," -H2O",sep="")
    blabel3<-paste(blabel1," -NH3",sep="")
    blabel4<-paste(blabel1," -98",sep="")
    
    ylabel1<-paste(c("y"),1:nby,sep="")
    ylabel2<-paste(ylabel1," -H2O",sep="")
    ylabel3<-paste(ylabel1," -NH3",sep="")
    ylabel4<-paste(ylabel1," -98",sep="")
    
    ####Obtain the decoy sequence fragment ions masses with the same b-y labels of  target peptide;
    fiby_label<-c(blabel1,blabel2,blabel3,blabel4,ylabel1,ylabel2,ylabel3,ylabel4)
    newbmass_razor<-newbmass[1:nby]
    newymass_razor<-newymass[1:nby]
    fibymass<-c(newbmass_razor,newbmass_razor-18.0152,newbmass_razor-17.0304,newbmass_razor-97.9815,
                newymass_razor,newymass_razor-18.0152,newymass_razor-17.0304,newymass_razor-97.9815)
    
    stime<-as.numeric(as.character(data_phsty$Start_Time[i_proba]))
    etime<-as.numeric(as.character(data_phsty$End_Time[i_proba]))
    newrt<-which(datapsm01$RT<=etime & datapsm01$RT>=stime)
    nnewrt<-length(newrt)
    char_sty<-gregexpr("S|T|Y",nomodseq[i_proba])
    nchar_sty<-length(char_sty[[1]])
    
    char_fi_mass<-as.character(data_phsty$FIMasses[i_proba])
    fimass_char<-strsplit(char_fi_mass,";")
    fimass<-as.numeric(fimass_char[[1]])
    nfimass<-length(fimass)
    fibymass_total<-c(fimass,fibymass)
    n_fibymass_total<-length(fibymass_total)
    
    fimass_ph1_intensity_proba<-vector()
    str_fi_intensity_proba<-vector()
    fimass_ph1_noise_proba<-vector()
    str_fi_noise_proba<-vector()
    for(sfi1_proba in 1:n_fibymass_total){
      for(nrt_proba in 1:nnewrt){
        nrt1_proba<-newrt[nrt_proba]
        char_mz_proba<-as.character(datapsm01$m.z[nrt1_proba])
        str_mz_proba<-strsplit(char_mz_proba,";")
        newstr_mz_proba<-as.numeric(str_mz_proba[[1]])
        char_intensity_proba<-as.character(datapsm01$Intensity[nrt1_proba])
        str_intensity_proba<-strsplit(char_intensity_proba,";")
        newstr_intensity_proba<-as.numeric(str_intensity_proba[[1]])
        char_noise_proba<-as.character(datapsm01$Noise[nrt1_proba])
        str_noise_proba<-strsplit(char_noise_proba,";")
        newstr_noise_proba<-as.numeric(str_noise_proba[[1]])
        
        sub_fion1_proba<-abs(newstr_mz_proba-fibymass_total[sfi1_proba])
        if(min(sub_fion1_proba)<=peakTolThreshold){
          sub_fion1_min_proba<-which.min(sub_fion1_proba)
          fimass_ph1_intensity_proba[nrt_proba]<-newstr_intensity_proba[sub_fion1_min_proba]
          fimass_ph1_noise_proba[nrt_proba]<-newstr_noise_proba[sub_fion1_min_proba]
        }else{
          fimass_ph1_intensity_proba[nrt_proba]<-0
          fimass_ph1_noise_proba[nrt_proba]<-0.0001
        }
      }
      newfimass_ph1_inten_proba<-paste(fimass_ph1_intensity_proba,collapse=";")
      str_fi_intensity_proba[sfi1_proba]<-newfimass_ph1_inten_proba
      newfimass_ph1_noise_proba<-paste(fimass_ph1_noise_proba,collapse=";")
      str_fi_noise_proba[sfi1_proba]<-newfimass_ph1_noise_proba
    }
    #get the XCorr
    spe_ref_XCorr_by_proba<-vector()
    found_frigmentions_by_intensity<-vector()
    found_frigmentions_by_intensity_log<-vector()
    ncor_by_proba<-1
    for(npby2 in (nfimass+1):n_fibymass_total){
      str_fi_intensitynew2_proba<-strsplit(str_fi_intensity_proba[npby2],";")
      str_fi_intensitynew3_proba<-as.numeric(str_fi_intensitynew2_proba[[1]])
      nnozero_str_fi_intensitynew3_1_proba<-sum(str_fi_intensitynew3_proba!=0)
      str_fi_noisenew2_proba<-strsplit(str_fi_noise_proba[npby2],";")
      str_fi_noisenew3_proba<-as.numeric(str_fi_noisenew2_proba[[1]])
      if(nnozero_str_fi_intensitynew3_1_proba<ppt){
        spe_ref_XCorr_by_proba[ncor_by_proba]<-0
        found_frigmentions_by_intensity[ncor_by_proba]<-0
        found_frigmentions_by_intensity_log[ncor_by_proba]<-0
      }else{
        cor_by1_proba<-vector()
        sum_str_fi_intensitynew1_proba<-vector()
        for(npby1_proba in 1:nfimass){
          str_fi_intensitynew_proba<-strsplit(str_fi_intensity_proba[npby1_proba],";")
          str_fi_intensitynew1_proba<-as.numeric(str_fi_intensitynew_proba[[1]])
          nnozero_str_fi_intensitynew1_proba<-sum(str_fi_intensitynew1_proba!=0)
          if(nnozero_str_fi_intensitynew1_proba<ppt){
            cor_by1_proba[npby1_proba]<-0
            sum_str_fi_intensitynew1_proba[npby1_proba]<-0
          }else{
            cor_by1_proba[npby1_proba]<-sum(str_fi_intensitynew1_proba*str_fi_intensitynew3_proba)/sqrt(sum(str_fi_intensitynew1_proba^2)*sum(str_fi_intensitynew3_proba^2))
            sum_str_fi_intensitynew1_proba[npby1_proba]<-log10(sum(str_fi_intensitynew1_proba))
          }
        }
        if(all(cor_by1_proba==0)){
          spe_ref_XCorr_by_proba[ncor_by_proba]<-0
          found_frigmentions_by_intensity[ncor_by_proba]<-0
          found_frigmentions_by_intensity_log[ncor_by_proba]<-0
        }else{
          cor_by_proba_weight<-sum(cor_by1_proba*sum_str_fi_intensitynew1_proba)/sum(sum_str_fi_intensitynew1_proba)
          cor_by_noise_proba<-str_fi_intensitynew3_proba/str_fi_noisenew3_proba
          ncor_by_noise_proba<-sum(cor_by_noise_proba>0)
          cor_by_noise_index_proba<-which(cor_by_noise_proba>0)
          cor_by_noise_sum_proba<-sum(cor_by_noise_proba[cor_by_noise_index_proba])
          spe_ref_XCorr_by_proba[ncor_by_proba]<-cor_by_proba_weight*(1-exp(-(cor_by_noise_sum_proba-ncor_by_noise_proba)))
          
          found_frigmentions_by_intensity_log[ncor_by_proba]<-log10(sum(str_fi_intensitynew3_proba))
          found_frigmentions_by_intensity[ncor_by_proba]<-sum(str_fi_intensitynew3_proba)
        }
      }
      ncor_by_proba<-ncor_by_proba+1
    }
    spe_ref_XCorr_total_sum[i_proba]<-sum(spe_ref_XCorr_by_proba)
    spe_ref_XCorr_total_sum_filter[i_proba]<-sum(spe_ref_XCorr_by_proba[which(spe_ref_XCorr_by_proba>=CorThreshold)])
    found_frigmentions_sum_intensity_log[i_proba]<-sum(found_frigmentions_by_intensity_log)
    found_frigmentions_sum_intensity_log_filter[i_proba]<-sum(found_frigmentions_by_intensity_log[which(spe_ref_XCorr_by_proba>=CorThreshold)])
    found_frigmentions_sum_intensity_total_filter[i_proba]<-sum(found_frigmentions_by_intensity[which(spe_ref_XCorr_by_proba>=CorThreshold)])
    Number_FragIons_Total[i_proba]<-sum(spe_ref_XCorr_by_proba>0)
    number_FragIons_aboveCorthreshold[i_proba]<-sum(spe_ref_XCorr_by_proba>=CorThreshold)
    found_frigmentions_sum_intensity_log_filter_average[i_proba]<-sum(found_frigmentions_by_intensity_log[which(spe_ref_XCorr_by_proba>=CorThreshold)])/sum(spe_ref_XCorr_by_proba>=CorThreshold)
    #WeightedCorr_total_average[i_proba]<-sum(spe_ref_XCorr_by_proba[which(spe_ref_XCorr_by_proba>=CorThreshold)])/sum(spe_ref_XCorr_by_proba>=CorThreshold)
    
    #setTxtProgressBar(pb2,i_proba)
  }
  data_phsty_after_proba<-cbind(data_phsty,Spe_Ref_XCorr_Total_Sum=spe_ref_XCorr_total_sum,
                                Spe_Ref_XCorr_Total_Sum_Filter=spe_ref_XCorr_total_sum_filter,
                                Found_Frigmentions_Sum_Intensity_Log=found_frigmentions_sum_intensity_log,
                                Found_Frigmentions_Sum_Intensity_Log_Filter=found_frigmentions_sum_intensity_log_filter,
                                Number_FragIons_Total=Number_FragIons_Total,
                                Number_FragIons_aboveCorthreshold=number_FragIons_aboveCorthreshold,
                                Found_Frigmentions_Sum_Intensity_Log_Filter_Average=found_frigmentions_sum_intensity_log_filter_average,
                                Found_Frigmentions_Sum_Intensity_Total_Filter=found_frigmentions_sum_intensity_total_filter)
  PSite_Change_FileName<-paste("Add_PSite_ChangedFiMasses_",NormalFileName,sep="")
  write.csv(data_phsty_after_proba,row.names=F,file=PSite_Change_FileName)
  #close(pb2)
}

#####05. Build the probable and decoy library

BDL<-function(ProbaTargetFileName,decoymethod="reverse",PSitechange=FALSE){
  data_phsty<-read.csv(ProbaTargetFileName,head=T,stringsAsFactors=F,sep=",")
  ndata_phsty<-length(data_phsty$File_Name)
  pro_decoy<-rep("Decoys",ndata_phsty)
  decoy_seq<-vector()
  modseq1<-as.character(data_phsty$Peptide_Modified_Sequence)
  modseq1<-gsub("\\[\\+42\\]","a",modseq1)
  modseq1<-gsub("\\[\\+80\\]","p",modseq1)
  modseq1<-gsub("\\[\\+16\\]","o",modseq1)
  modseq1<-gsub("\\[\\+57\\]","c",modseq1)
  modseq1<-gsub("\\[\\+10\\]","r",modseq1)
  modseq1<-gsub("\\[\\+8\\]","l",modseq1)
  pb<- txtProgressBar(min=1,max=ndata_phsty,char="=>")
  for(idecoy in 1:ndata_phsty){
    newmodseq1<-strsplit(modseq1[idecoy],"")
    nnewmodseq1<-length(newmodseq1[[1]])
    
    newmodseq1_rev<-vector()
    for(i in 1:nnewmodseq1){
      if(newmodseq1[[1]][i]=="a" | newmodseq1[[1]][i]=="p" | newmodseq1[[1]][i]=="o"| newmodseq1[[1]][i]=="c"|newmodseq1[[1]][i]=="r"|newmodseq1[[1]][i]=="l"){
        newmodseq1_rev[i-1]<-paste(newmodseq1[[1]][i-1],newmodseq1[[1]][i],sep="")
      }else{
        newmodseq1_rev[i]<-newmodseq1[[1]][i]
      }
    }
    newmodseq2_rev<-na.omit(newmodseq1_rev)
    nnewmodseq2_rev<-length(newmodseq2_rev)
    if(PSitechange==TRUE){
      char_rev_ph_index<-0
    }else{
      char_rev_ph<-gregexpr("p",newmodseq2_rev[1:(nnewmodseq2_rev-1)])
      char_rev_ph_index<-which(charrevph>0)
    }
    
    #newmodseq2_rev[char_rev_ph_index]<-NA
    #newmodseq3_rev<-na.omit(newmodseq2_rev)
    
    newmodseq3_rev<-newmodseq2_rev[1:(nnewmodseq2_rev-1)]
    newmodseq3_rev[char_rev_ph_index]<-NA
    newmodseq3_rev<-na.omit(newmodseq3_rev)
    if(decoymethod=="reverse"){
      decoymethod1<-"reverse_"
      newmodseq4_rev<-rev(newmodseq3_rev)
    }
    if(decoymethod=="shuffle"){
      decoymethod1<-"shuffle_"
      newmodseq4_rev<-sample(newmodseq3_rev)
    }
    #newmodseq2_rev_decoy<-vector(length=nnewmodseq2_rev)
    
    newmodseq2_rev_decoy<-vector(length=(nnewmodseq2_rev-1))
    newmodseq2_rev_decoy[char_rev_ph_index]<-newmodseq2_rev[char_rev_ph_index]
    newmodseq2_rev_decoy_index<-which(newmodseq2_rev_decoy=="FALSE")
    newmodseq2_rev_decoy[newmodseq2_rev_decoy_index]<-newmodseq4_rev
    
    newrev_seq<-paste(newmodseq2_rev_decoy,collapse="")
    newrev_seq_intact<-paste(newrev_seq,newmodseq2_rev[nnewmodseq2_rev],sep="")
    newrev_seq1<-gsub("a","[+42]",newrev_seq_intact)
    newrev_seq1<-gsub("p","[+80]",newrev_seq1)
    newrev_seq1<-gsub("o","[+16]",newrev_seq1)
    newrev_seq1<-gsub("c","[+57]",newrev_seq1)
    newrev_seq1<-gsub("r","[+10]",newrev_seq1)
    newrev_seq1<-gsub("l","[+8]",newrev_seq1)
    decoy_seq[idecoy]<-newrev_seq1
    
    setTxtProgressBar(pb,idecoy)
  }
  
  newdata_phsty<-data.frame(File_Name=rep(as.character(data_phsty$File_Name),2),
                            Protein_Name=c(as.character(data_phsty$Protein_Name),pro_decoy),
                            Peptide_Modified_Sequence=c(as.character(data_phsty$Peptide_Modified_Sequence),decoy_seq),
                            Precursor_Mass=rep(data_phsty$Precursor_Mass,2),
                            Precursor_Charge=rep(data_phsty$Precursor_Charge,2),
                            Fragment_Ions=rep(data_phsty$Fragment_Ions,2),
                            FIMasses=rep(data_phsty$FIMasses,2),
                            Rentention_Time=rep(data_phsty$Rentention_Time,2),
                            Start_Time=rep(data_phsty$Start_Time,2),
                            End_Time=rep(data_phsty$End_Time,2),
                            Serial_Number_Extension=rep(data_phsty$Serial_Number_Extension,2),
                            Spe_Ref_XCorr_Total_Sum=c(data_phsty$Spe_Ref_XCorr_Total_Sum,rep(-1,ndata_phsty)),
                            Spe_Ref_XCorr_Total_Sum_Filter=c(data_phsty$Spe_Ref_XCorr_Total_Sum_Filter,rep(-1,ndata_phsty)),
                            Found_Frigmentions_Sum_Intensity_Log=c(data_phsty$Found_Frigmentions_Sum_Intensity_Log,rep(-1,ndata_phsty)),
                            Found_Frigmentions_Sum_Intensity_Log_Filter=c(data_phsty$Found_Frigmentions_Sum_Intensity_Log_Filter,rep(-1,ndata_phsty)),
                            Number_FragIons_Total=c(data_phsty$Number_FragIons_Total,rep(-1,ndata_phsty)),
                            Number_FragIons_aboveCorthreshold_Total=c(data_phsty$Number_FragIons_aboveCorthreshold,rep(-1,ndata_phsty)),
                            Found_Frigmentions_Sum_Intensity_Log_Filter_Average=c(data_phsty$Found_Frigmentions_Sum_Intensity_Log_Filter_Average,rep(-1,ndata_phsty)),
                            Found_Frigmentions_Sum_Intensity_Total_Filter=c(data_phsty$Found_Frigmentions_Sum_Intensity_Total_Filter,rep(-1,ndata_phsty)))
  if(PSitechange==TRUE){
    targetfn1<-paste("TargetDecoy_",decoymethod1,sep="")
    decoyfn1<-paste("Decoy_",decoymethod1,sep="")
  }else{
    targetfn1<-paste("TargetDecoy_phositenc_",decoymethod1,sep="")
    decoyfn1<-paste("Decoy_phositenc_",decoymethod1,sep="")
  }
  TDFileName<-paste(targetfn1,ProbaTargetFileName,sep="") #n(ph) != n(STY)
  write.csv(newdata_phsty,row.names=F,file=TDFileName)
  
  decoyfn<-paste(decoyfn1,ProbaTargetFileName,sep="")
  data_phsty_decoy<-newdata_phsty[(ndata_phsty+1):(2*ndata_phsty),]
  write.csv(data_phsty_decoy,row.names=F,file=decoyfn)
  
  close(pb)
}

########06. Evaluate the phosphosites

EDIA<-function(TDFileName,FinalFileName,RawFileName,IsoWindowFileName,CorThreshold,
               ppt=3,peakTolThreshold=0.05){
  data_rawms2<-read.table(RawFileName,head=T,stringsAsFactors=F,sep="\t")
  nrawms2<-length(data_rawms2$ScanNum)
  data_isolationw<-read.csv(IsoWindowFileName,head=T,stringsAsFactors=F,sep=",")
  ndata_isolationw<-length(data_isolationw[[1]])
  niw_times<-nrawms2%/%ndata_isolationw
  newdata_isolationw<-rep(data_isolationw$Center,niw_times)
  data_isolationw_after<-data.frame(IW_Center=newdata_isolationw) #Isolation Window Center
  n_isow<-niw_times*ndata_isolationw
  newdata_rawms2<-cbind(data_rawms2[1:n_isow,],data_isolationw_after)
  
  modifiedaa<-hash(Hyd=1.0079,Car=12.0107,Nit=14.0067,Oxi=15.9994,Pho=30.9738,Sul=32.065,
                   ac=42.01056,ox=15.99491,ph=79.9663,ca=57.0215,ly=8.0142,ar=10.00827,
                   A=71.03711,R=156.10111,N=114.04293,D=115.02694,C=103.00919,
                   E=129.04259,Q=128.05858,G=57.02146,H=137.05891,I=113.08406,
                   L=113.08406,K=128.09496,M=131.04049,F=147.06841,P=97.05276,
                   S=87.03203,T=101.04768,W=186.07931,Y=163.06333,V=99.06841)
  
  data_phsty<-read.csv(TDFileName,head=T,stringsAsFactors=F,sep=",")
  nms<-length(data_phsty[[1]])/2
  modseq<-data_phsty$Peptide_Modified_Sequence
  modseq<-gsub("\\[\\+42\\]","a",modseq)
  modseq<-gsub("\\[\\+80\\]","p",modseq)
  modseq<-gsub("\\[\\+16\\]","o",modseq)
  modseq<-gsub("\\[\\+57\\]","c",modseq)
  modseq<-gsub("\\[\\+10\\]","r",modseq)
  modseq<-gsub("\\[\\+8\\]","l",modseq)
  modseq_ac<-gsub("p|o|c|r|l","",modseq)
  modseq_ph<-gsub("a|o|c|r|l","",modseq)
  modseq_ox<-gsub("a|p|c|r|l","",modseq)
  modseq_ca<-gsub("a|p|o|r|l","",modseq)
  modseq_ar<-gsub("a|p|o|l|c","",modseq)
  modseq_ly<-gsub("a|p|o|c|r","",modseq)
  nomodseq<-gsub("a|o|c|p|r|l","",modseq)
  
  spe_ref_XCorr_sum<-vector()
  spe_ACorr_sum<-vector()
  ref_ACorr_sum<-vector()
  spe_fi_num<-vector()
  spe_fi_intensity_sum<-vector()
  ref_fi_num<-vector()
  ref_fi_intensity_sum<-vector()
  max_cor_by<-vector()
  pb<-txtProgressBar(min=1,max=nms,char="=>")
  for(i in 1:nms){
    iw_compare1<-abs(data_phsty$Precursor_Mass[i]-data_isolationw$Center[1:(ndata_isolationw/2)])
    iw_compare_index1<-which.min(iw_compare1)
    sub_ms2_iw_center_precursor1<-data_isolationw$Center[1:(ndata_isolationw/2)][iw_compare_index1]
    iw_compare2<-abs(data_phsty$Precursor_Mass[i]-data_isolationw$Center[(ndata_isolationw/2+1):ndata_isolationw])
    iw_compare_index2<-which.min(iw_compare2)
    sub_ms2_iw_center_precursor2<-data_isolationw$Center[(ndata_isolationw/2+1):ndata_isolationw][iw_compare_index2]
    
    iw_each_index<-vector()
    iw_each_index1<-which(newdata_rawms2$IW_Center==sub_ms2_iw_center_precursor1)
    iw_each_index2<-which(newdata_rawms2$IW_Center==sub_ms2_iw_center_precursor2)
    jw<-1
    for(iw in 1:length(iw_each_index1)){
      iw_each_index[jw]<-iw_each_index1[iw]
      iw_each_index[jw+1]<-iw_each_index2[iw]
      jw<-jw+2
    }
    iw_each<-newdata_rawms2[iw_each_index,]
    datapsm01<-iw_each
    #for(jisow in 1:ndata_isolationw){
    #  if(data_phsty$Precursor_Mass[i]<data_isolationw$End[jisow] &data_phsty$Precursor_Mass[i]>=data_isolationw$Start[jisow]){
    #    sub_ms2_iw_center_precursor<-data_isolationw$Center[jisow]
    #    break
    #  }
    #}
    #iw_each_index<-which(newdata_rawms2$IW_Center==sub_ms2_iw_center_precursor)
    #iw_each<-newdata_rawms2[iw_each_index,]
    #datapsm01<-iw_each
    ####
    ####b,ymass decoy
    bmass_decoy<-vector()
    ymass_decoy<-vector()
    noms_split<-strsplit(nomodseq[i+nms],"")
    nnoms_split<-length(noms_split[[1]])
    bmass_decoy[1]<-modifiedaa[[noms_split[[1]][1]]]+modifiedaa$Hyd
    bb<-bmass_decoy[1]
    ymass_decoy[1]<-modifiedaa[[noms_split[[1]][nnoms_split]]]+modifiedaa$Hyd*3+modifiedaa$Oxi
    yy<-ymass_decoy[1]
    for(jb in 2:nnoms_split){
      bb<-bb+modifiedaa[[noms_split[[1]][jb]]]
      bmass_decoy[jb]<-bb
    }
    for(jy in 2:nnoms_split){
      yy<-yy+modifiedaa[[noms_split[[1]][nnoms_split-jy+1]]]
      ymass_decoy[jy]<-yy
    }
    
    newbmass_decoy<-bmass_decoy
    newymass_decoy<-ymass_decoy
    nby<-nnoms_split
    char_modseq_ac<-gregexpr("a",modseq_ac[i+nms])
    char_modseq_ox<-gregexpr("o",modseq_ox[i+nms])
    char_modseq_ph<-gregexpr("p",modseq_ph[i+nms])
    char_modseq_ca<-gregexpr("c",modseq_ca[i+nms])
    char_modseq_ar<-gregexpr("r",modseq_ar[i+nms])
    char_modseq_ly<-gregexpr("l",modseq_ly[i+nms])
    nseq_ac<-sum(char_modseq_ac[[1]]>=1)
    nseq_ox<-sum(char_modseq_ox[[1]]>=1)
    nseq_ph<-sum(char_modseq_ph[[1]]>=1)
    nseq_ca<-sum(char_modseq_ca[[1]]>=1)
    nseq_ar<-sum(char_modseq_ar[[1]]>=1)
    nseq_ly<-sum(char_modseq_ly[[1]]>=1)
    if(nseq_ac>=1){
      for(acj in 1:nseq_ac){
        nacj<-char_modseq_ac[[1]][acj]-acj
        newbmass_decoy[nacj:nby]<-newbmass_decoy[nacj:nby]+modifiedaa$ac
        newymass_decoy[(nby-nacj+1):nby]<-newymass_decoy[(nby-nacj+1):nby]+modifiedaa$ac
      }
    }
    if(nseq_ox>=1){
      for(oxj in 1:nseq_ox){
        noxj<-char_modseq_ox[[1]][oxj]-oxj
        newbmass_decoy[noxj:nby]<-newbmass_decoy[noxj:nby]+modifiedaa$ox
        newymass_decoy[(nby-noxj+1):nby]<-newymass_decoy[(nby-noxj+1):nby]+modifiedaa$ox
      }
    }
    if(nseq_ca>=1){
      for(caj in 1:nseq_ca){
        ncaj<-char_modseq_ca[[1]][caj]-caj
        newbmass_decoy[ncaj:nby]<-newbmass_decoy[ncaj:nby]+modifiedaa$ca
        newymass_decoy[(nby-ncaj+1):nby]<-newymass_decoy[(nby-ncaj+1):nby]+modifiedaa$ca
      }
    }
    if(nseq_ph>=1){
      for(phj in 1:nseq_ph){
        nphj<-char_modseq_ph[[1]][phj]-phj
        newbmass_decoy[nphj:nby]<-newbmass_decoy[nphj:nby]+modifiedaa$ph
        newymass_decoy[(nby-nphj+1):nby]<-newymass_decoy[(nby-nphj+1):nby]+modifiedaa$ph
      }
    }
    if(nseq_ar>=1){
      for(arj in 1:nseq_ar){
        narj<-char_modseq_ar[[1]][arj]-arj
        newbmass_decoy[narj:nby]<-newbmass_decoy[narj:nby]+modifiedaa$ar
        newymass_decoy[(nby-narj+1):nby]<-newymass_decoy[(nby-narj+1):nby]+modifiedaa$ar
      }
    }
    if(nseq_ly>=1){
      for(lyj in 1:nseq_ly){
        nlyj<-char_modseq_ly[[1]][lyj]-lyj
        newbmass_decoy[nlyj:nby]<-newbmass_decoy[nlyj:nby]+modifiedaa$ly
        newymass_decoy[(nby-nlyj+1):nby]<-newymass_decoy[(nby-nlyj+1):nby]+modifiedaa$ly
      }
    }
    newymass_decoy<-newymass_decoy-0.005
    yyrev_decoy<-rev(newymass_decoy)
    
    ####
    bmass<-vector()
    ymass<-vector()
    noms_split<-strsplit(nomodseq[i],"")
    nnoms_split<-length(noms_split[[1]])
    bmass[1]<-modifiedaa[[noms_split[[1]][1]]]+modifiedaa$Hyd
    bb<-bmass[1]
    ymass[1]<-modifiedaa[[noms_split[[1]][nnoms_split]]]+modifiedaa$Hyd*3+modifiedaa$Oxi
    yy<-ymass[1]
    for(jb in 2:nnoms_split){
      bb<-bb+modifiedaa[[noms_split[[1]][jb]]]
      bmass[jb]<-bb
    }
    for(jy in 2:nnoms_split){
      yy<-yy+modifiedaa[[noms_split[[1]][nnoms_split-jy+1]]]
      ymass[jy]<-yy
    }
    
    newbmass<-bmass
    newymass<-ymass
    nby<-nnoms_split
    char_modseq_ac<-gregexpr("a",modseq_ac[i])
    char_modseq_ox<-gregexpr("o",modseq_ox[i])
    char_modseq_ph<-gregexpr("p",modseq_ph[i])
    char_modseq_ca<-gregexpr("c",modseq_ca[i])
    char_modseq_ar<-gregexpr("r",modseq_ar[i])
    char_modseq_ly<-gregexpr("l",modseq_ly[i])
    nseq_ac<-sum(char_modseq_ac[[1]]>=1)
    nseq_ox<-sum(char_modseq_ox[[1]]>=1)
    nseq_ph<-sum(char_modseq_ph[[1]]>=1)
    nseq_ca<-sum(char_modseq_ca[[1]]>=1)
    nseq_ar<-sum(char_modseq_ar[[1]]>=1)
    nseq_ly<-sum(char_modseq_ly[[1]]>=1)
    if(nseq_ac>=1){
      for(acj in 1:nseq_ac){
        nacj<-char_modseq_ac[[1]][acj]-acj
        newbmass[nacj:nby]<-newbmass[nacj:nby]+modifiedaa$ac
        newymass[(nby-nacj+1):nby]<-newymass[(nby-nacj+1):nby]+modifiedaa$ac
      }
    }
    if(nseq_ox>=1){
      for(oxj in 1:nseq_ox){
        noxj<-char_modseq_ox[[1]][oxj]-oxj
        newbmass[noxj:nby]<-newbmass[noxj:nby]+modifiedaa$ox
        newymass[(nby-noxj+1):nby]<-newymass[(nby-noxj+1):nby]+modifiedaa$ox
      }
    }
    if(nseq_ca>=1){
      for(caj in 1:nseq_ca){
        ncaj<-char_modseq_ca[[1]][caj]-caj
        newbmass[ncaj:nby]<-newbmass[ncaj:nby]+modifiedaa$ca
        newymass[(nby-ncaj+1):nby]<-newymass[(nby-ncaj+1):nby]+modifiedaa$ca
      }
    }
    if(nseq_ph>=1){
      for(phj in 1:nseq_ph){
        nphj<-char_modseq_ph[[1]][phj]-phj
        newbmass[nphj:nby]<-newbmass[nphj:nby]+modifiedaa$ph
        newymass[(nby-nphj+1):nby]<-newymass[(nby-nphj+1):nby]+modifiedaa$ph
      }
    }
    if(nseq_ar>=1){
      for(arj in 1:nseq_ar){
        narj<-char_modseq_ar[[1]][arj]-arj
        newbmass[narj:nby]<-newbmass[narj:nby]+modifiedaa$ar
        newymass[(nby-narj+1):nby]<-newymass[(nby-narj+1):nby]+modifiedaa$ar
      }
    }
    if(nseq_ly>=1){
      for(lyj in 1:nseq_ly){
        nlyj<-char_modseq_ly[[1]][lyj]-lyj
        newbmass[nlyj:nby]<-newbmass[nlyj:nby]+modifiedaa$ly
        newymass[(nby-nlyj+1):nby]<-newymass[(nby-nlyj+1):nby]+modifiedaa$ly
      }
    }
    newymass<-newymass-0.005
    yyrev<-rev(newymass)
    
    ####
    blabel1<-paste(c("b"),1:nby,sep="")
    blabel2<-paste(blabel1," -H2O",sep="")
    blabel3<-paste(blabel1," -NH3",sep="")
    blabel4<-paste(blabel1," -98",sep="")
    
    ylabel1<-paste(c("y"),1:nby,sep="")
    ylabel2<-paste(ylabel1," -H2O",sep="")
    ylabel3<-paste(ylabel1," -NH3",sep="")
    ylabel4<-paste(ylabel1," -98",sep="")
    
    ####Obtain the decoy sequence fragment ions masses with the same b-y labels of  target peptide;
    fiby_label<-c(blabel1,blabel4,ylabel1,ylabel4)
    fibymass_decoy<-hash(fiby_label,c(newbmass_decoy,newbmass_decoy-97.9815,
                                      newymass_decoy,newymass_decoy-97.9815))
    char_fi<-as.character(data_phsty$Fragment_Ions[i])
    fi_char<-strsplit(char_fi,";")
    nfi_char<-length(fi_char[[1]])
    fimass_decoy<-vector()
    fi_not_match<-vector()
    jfichar<-1
    jfi_char_not_match<-1
    for(ifichar in 1:nfi_char){
      fi_char_match_index<-which(fiby_label==fi_char[[1]][ifichar])
      nfi_char_match<-sum(fi_char_match_index)
      if(nfi_char_match>0){
        fimass_decoy[jfichar]<-fibymass_decoy[[fi_char[[1]][ifichar]]]
        jfichar<-jfichar+1
      }else{
        fi_not_match[jfi_char_not_match]<-ifichar
        jfi_char_not_match<-jfi_char_not_match+1
      }
    }
    #### If the decoy peptide does not contain some b-y ions, remove also those from target peptide;
    char_fi_mass<-as.character(data_phsty$FIMasses[i])
    fimass_char<-strsplit(char_fi_mass,";")
    fimass1<-as.numeric(fimass_char[[1]])
    nfi_not_match<-length(fi_not_match)
    if(nfi_not_match>0){
      fimass<-fimass1[-fi_not_match]
    }else{
      fimass<-fimass1
    }
    nfimass<-length(fimass)
    
    stime<-as.numeric(as.character(data_phsty$Start_Time[i]))
    etime<-as.numeric(as.character(data_phsty$End_Time[i]))
    newrt<-which(datapsm01$RT<=etime & datapsm01$RT>=stime)
    nnewrt<-length(newrt)
    char_sty<-gregexpr("S|T|Y",nomodseq[i])
    nchar_sty<-length(char_sty[[1]])
    ####
    spe_ref_XCorr_sum_ph<-vector()
    spe_fi_num_ph<-vector()
    spe_fi_intensity_sum_ph<-vector()
    max_cor_by_ph<-vector()
    
    spe_ref_XCorr_sum_ph_decoy<-vector()
    spe_fi_num_ph_decoy<-vector()
    spe_fi_intensity_sum_ph_decoy<-vector()
    max_cor_by_ph_decoy<-vector()
    
    #####
    for(phk in 1:nseq_ph){
      nphk<-char_modseq_ph[[1]][phk]-phk
      
      if(nphk==char_sty[[1]][1]){
        phk_bmass_after11_start<-newbmass[nphk:(char_sty[[1]][2]-1)]
        phk_bmass_after12_start<-phk_bmass_after11_start-18.0152
        phk_bmass_after13_start<-phk_bmass_after11_start-17.0304
        phk_bmass_after1_start<-c(phk_bmass_after11_start,phk_bmass_after12_start,
                                  phk_bmass_after13_start)
        phk_ymass_ahead11_start<-yyrev[(nphk+1):char_sty[[1]][2]]
        phk_ymass_ahead12_start<-phk_ymass_ahead11_start-18.0152
        phk_ymass_ahead13_start<-phk_ymass_ahead11_start-17.0304
        phk_ymass_ahead1_start<-c(phk_ymass_ahead11_start,phk_ymass_ahead12_start,
                                  phk_ymass_ahead13_start)
        fi_bfmass_start<-c(fimass,phk_bmass_after1_start)
        fi_yhmass_start<-c(fimass,phk_ymass_ahead1_start)
        nfi_bfmass_start<-length(fi_bfmass_start)
        nfi_yhmass_start<-length(fi_yhmass_start)
        
        fimass_ph1_intensity_start<-vector()
        str_fi_intensity_start<-vector()
        fimass_ph1_noise_start<-vector()
        str_fi_noise_start<-vector()
        
        ####
        phk_bmass_after_decoy11_start<-newbmass_decoy[nphk:(char_sty[[1]][2]-1)]
        phk_bmass_after_decoy12_start<-phk_bmass_after_decoy11_start-18.0152
        phk_bmass_after_decoy13_start<-phk_bmass_after_decoy11_start-17.0304
        phk_bmass_after_decoy1_start<-c(phk_bmass_after_decoy11_start,phk_bmass_after_decoy12_start,
                                        phk_bmass_after_decoy13_start)
        phk_ymass_ahead_decoy11_start<-yyrev_decoy[(nphk+1):char_sty[[1]][2]]
        phk_ymass_ahead_decoy12_start<-phk_ymass_ahead_decoy11_start-18.0152
        phk_ymass_ahead_decoy13_start<-phk_ymass_ahead_decoy11_start-17.0304
        phk_ymass_ahead_decoy1_start<-c(phk_ymass_ahead_decoy11_start,phk_ymass_ahead_decoy12_start,
                                        phk_ymass_ahead_decoy13_start)
        fi_bfmass_decoy_start<-c(fimass_decoy,phk_bmass_after_decoy1_start)
        fi_yhmass_decoy_start<-c(fimass_decoy,phk_ymass_ahead_decoy1_start)
        nfi_bfmass_decoy_start<-length(fi_bfmass_decoy_start)
        nfi_yhmass_decoy_start<-length(fi_yhmass_decoy_start)
        
        fimass_ph1_intensity_decoy_start<-vector()
        str_fi_intensity_decoy_start<-vector()
        fimass_ph1_noise_decoy_start<-vector()
        str_fi_noise_decoy_start<-vector()
        
        ####
        for(sfi1_start in 1:nfi_bfmass_start){
          for(nrt_start in 1:nnewrt){
            nrt1_start<-newrt[nrt_start]
            char_mz_start<-as.character(datapsm01$m.z[nrt1_start])
            str_mz_start<-strsplit(char_mz_start,";")
            newstr_mz_start<-as.numeric(str_mz_start[[1]])
            char_intensity_start<-as.character(datapsm01$Intensity[nrt1_start])
            str_intensity_start<-strsplit(char_intensity_start,";")
            newstr_intensity_start<-as.numeric(str_intensity_start[[1]])
            char_noise_start<-as.character(datapsm01$Noise[nrt1_start])
            str_noise_start<-strsplit(char_noise_start,";")
            newstr_noise_start<-as.numeric(str_noise_start[[1]])
            
            sub_fion1_start<-abs(newstr_mz_start-fi_bfmass_start[sfi1_start])
            if(min(sub_fion1_start)<=peakTolThreshold){
              sub_fion1_min_start<-which.min(sub_fion1_start)
              fimass_ph1_intensity_start[nrt_start]<-newstr_intensity_start[sub_fion1_min_start]
              fimass_ph1_noise_start[nrt_start]<-newstr_noise_start[sub_fion1_min_start]
            }else{
              fimass_ph1_intensity_start[nrt_start]<-0
              fimass_ph1_noise_start[nrt_start]<-0.0001
            }
            ####
            sub_fion1_decoy_start<-abs(newstr_mz_start-fi_bfmass_decoy_start[sfi1_start])
            if(min(sub_fion1_decoy_start)<=peakTolThreshold){
              sub_fion1_min_decoy_start<-which.min(sub_fion1_decoy_start)
              fimass_ph1_intensity_decoy_start[nrt_start]<-newstr_intensity_start[sub_fion1_min_decoy_start]
              fimass_ph1_noise_decoy_start[nrt_start]<-newstr_noise_start[sub_fion1_min_decoy_start]
            }else{
              fimass_ph1_intensity_decoy_start[nrt_start]<-0
              fimass_ph1_noise_decoy_start[nrt_start]<-0.0001
            }
          }
          newfimass_ph1_inten_start<-paste(fimass_ph1_intensity_start,collapse=";")
          str_fi_intensity_start[sfi1_start]<-newfimass_ph1_inten_start
          newfimass_ph1_noise_start<-paste(fimass_ph1_noise_start,collapse=";")
          str_fi_noise_start[sfi1_start]<-newfimass_ph1_noise_start
          ####
          newfimass_ph1_inten_decoy_start<-paste(fimass_ph1_intensity_decoy_start,collapse=";")
          str_fi_intensity_decoy_start[sfi1_start]<-newfimass_ph1_inten_decoy_start
          newfimass_ph1_noise_decoy_start<-paste(fimass_ph1_noise_decoy_start,collapse=";")
          str_fi_noise_decoy_start[sfi1_start]<-newfimass_ph1_noise_decoy_start
        }
        
        #####Target_b ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bf_start<-vector()
        spe_bf_intensity_start<-vector()
        ncor_bf_start<-1
        for(npbf2_start in (nfimass+1):nfi_bfmass_start){
          str_fi_intensitynew2_start<-strsplit(str_fi_intensity_start[npbf2_start],";")
          str_fi_intensitynew3_start<-as.numeric(str_fi_intensitynew2_start[[1]])
          nnozero_str_fi_intensitynew3_1_start<-sum(str_fi_intensitynew3_start!=0)
          str_fi_noisenew2_start<-strsplit(str_fi_noise_start[npbf2_start],";")
          str_fi_noisenew3_start<-as.numeric(str_fi_noisenew2_start[[1]])
          if(nnozero_str_fi_intensitynew3_1_start<ppt){
            spe_ref_XCorr_bf_start[ncor_bf_start]<-0
            spe_bf_intensity_start[ncor_bf_start]<-0
          }else{
            cor_bf1_start<-vector()
            sum_str_fi_intensitynew1_start<-vector()
            for(npbf1_start in 1:nfimass){
              str_fi_intensitynew_start<-strsplit(str_fi_intensity_start[npbf1_start],";")
              str_fi_intensitynew1_start<-as.numeric(str_fi_intensitynew_start[[1]])
              nnozero_str_fi_intensitynew1_start<-sum(str_fi_intensitynew1_start!=0)
              if(nnozero_str_fi_intensitynew1_start<ppt){
                cor_bf1_start[npbf1_start]<-0
                sum_str_fi_intensitynew1_start[npbf1_start]<-0
              }else{
                cor_bf1_start[npbf1_start]<-sum(str_fi_intensitynew1_start*str_fi_intensitynew3_start)/sqrt(sum(str_fi_intensitynew1_start^2)*sum(str_fi_intensitynew3_start^2))
                sum_str_fi_intensitynew1_start[npbf1_start]<-log2(sum(str_fi_intensitynew1_start))
              }
            }
            if(all(cor_bf1_start==0)){
              spe_ref_XCorr_bf_start[ncor_bf_start]<-0
              spe_bf_intensity_start[ncor_bf_start]<-0
            }else{
              cor_bf_start_weight<-sum(cor_bf1_start)#sum(cor_bf1_start*sum_str_fi_intensitynew1_start)/sum(sum_str_fi_intensitynew1_start)
              #cor_bf_noise_start<-str_fi_intensitynew3_start/str_fi_noisenew3_start
              #ncor_bf_noise_start<-sum(cor_bf_noise_start>0)
              #cor_bf_noise_index_start<-which(cor_bf_noise_start>0)
              #cor_bf_noise_sum_start<-sum(cor_bf_noise_start[cor_bf_noise_index_start])
              if(cor_bf_start_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bf_start[ncor_bf_start]<-cor_bf_start_weight
                nozero_str_fi_intensitynew3_start<-log2(sum(str_fi_intensitynew3_start))
                spe_bf_intensity_start[ncor_bf_start]<-nozero_str_fi_intensitynew3_start
              }else{
                spe_ref_XCorr_bf_start[ncor_bf_start]<-0
                spe_bf_intensity_start[ncor_bf_start]<-0
              }
              #*(1-exp(-(cor_bf_noise_sum_start-ncor_bf_noise_start)))
            }
          }
          ncor_bf_start<-ncor_bf_start+1
        }
        spe_ref_XCorr_bf_sum_start<-sum(spe_ref_XCorr_bf_start)
        spe_bf_intensity_sum_start<-sum(spe_bf_intensity_start)
        
        ####
        spe_bf_num_start<-sum(spe_ref_XCorr_bf_start!=0)
        
        if(max(spe_ref_XCorr_bf_start)>0){
          max_cor_index_start<-which.max(spe_ref_XCorr_bf_start)
          if(max_cor_index_start<=(nfi_bfmass_start-nfimass)/3){
            max_cor_by_ph_bf_start<-blabel1[nphk+max_cor_index_start-1]
          }
          else if(max_cor_index_start<=2*(nfi_bfmass_start-nfimass)/3 & max_cor_index_start>(nfi_bfmass_start-nfimass)/3){
            max_cor_index1_start<-max_cor_index_start-(nfi_bfmass_start-nfimass)/3
            max_cor_by_ph_bf_start<-blabel2[nphk+max_cor_index1_start-1]
          }
          else{
            max_cor_index2_start<-max_cor_index_start-2*(nfi_bfmass_start-nfimass)/3
            max_cor_by_ph_bf_start<-blabel3[nphk+max_cor_index2_start-1]
          }
        }else{
          max_cor_by_ph_bf_start<-"NA"
        }
        
        #####Decoy_b ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bf_decoy_start<-vector()
        spe_bf_intensity_decoy_start<-vector()
        ncor_bf_decoy_start<-1
        for(npbf2_start in (nfimass+1):nfi_bfmass_start){
          str_fi_intensitynew2_start<-strsplit(str_fi_intensity_decoy_start[npbf2_start],";")
          str_fi_intensitynew3_start<-as.numeric(str_fi_intensitynew2_start[[1]])
          nnozero_str_fi_intensitynew3_1_start<-sum(str_fi_intensitynew3_start!=0)
          str_fi_noisenew2_start<-strsplit(str_fi_noise_decoy_start[npbf2_start],";")
          str_fi_noisenew3_start<-as.numeric(str_fi_noisenew2_start[[1]])
          if(nnozero_str_fi_intensitynew3_1_start<ppt){
            spe_ref_XCorr_bf_decoy_start[ncor_bf_decoy_start]<-0
            spe_bf_intensity_decoy_start[ncor_bf_decoy_start]<-0
          }else{
            cor_bf1_start<-vector()
            sum_str_fi_intensitynew1_start<-vector()
            for(npbf1_start in 1:nfimass){
              str_fi_intensitynew_start<-strsplit(str_fi_intensity_decoy_start[npbf1_start],";")
              str_fi_intensitynew1_start<-as.numeric(str_fi_intensitynew_start[[1]])
              nnozero_str_fi_intensitynew1_start<-sum(str_fi_intensitynew1_start!=0)
              if(nnozero_str_fi_intensitynew1_start<ppt){
                cor_bf1_start[npbf1_start]<-0
                sum_str_fi_intensitynew1_start[npbf1_start]<-0
              }else{
                cor_bf1_start[npbf1_start]<-sum(str_fi_intensitynew1_start*str_fi_intensitynew3_start)/sqrt(sum(str_fi_intensitynew1_start^2)*sum(str_fi_intensitynew3_start^2))
                sum_str_fi_intensitynew1_start[npbf1_start]<-log2(sum(str_fi_intensitynew1_start))
              }
            }
            if(all(cor_bf1_start==0)){
              spe_ref_XCorr_bf_decoy_start[ncor_bf_decoy_start]<-0
              spe_bf_intensity_decoy_start[ncor_bf_decoy_start]<-0
            }else{
              cor_bf_start_weight<-sum(cor_bf1_start)#*sum_str_fi_intensitynew1_start)/sum(sum_str_fi_intensitynew1_start)
              #cor_bf_noise_start<-str_fi_intensitynew3_start/str_fi_noisenew3_start
              #ncor_bf_noise_start<-sum(cor_bf_noise_start>0)
              #cor_bf_noise_index_start<-which(cor_bf_noise_start>0)
              #cor_bf_noise_sum_start<-sum(cor_bf_noise_start[cor_bf_noise_index_start])
              if(cor_bf_start_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bf_decoy_start[ncor_bf_decoy_start]<-cor_bf_start_weight
                nozero_str_fi_intensitynew3_start<-log2(sum(str_fi_intensitynew3_start))#*nnozero_str_fi_intensitynew3_start
                spe_bf_intensity_decoy_start[ncor_bf_decoy_start]<-nozero_str_fi_intensitynew3_start
              }else{
                spe_ref_XCorr_bf_decoy_start[ncor_bf_decoy_start]<-0
                spe_bf_intensity_decoy_start[ncor_bf_decoy_start]<-0
              }
              #*(1-exp(-(cor_bf_noise_sum_start-ncor_bf_noise_start)))
            }
          }
          ncor_bf_decoy_start<-ncor_bf_decoy_start+1
        }
        spe_ref_XCorr_bf_sum_decoy_start<-sum(spe_ref_XCorr_bf_decoy_start)
        spe_bf_intensity_sum_decoy_start<-sum(spe_bf_intensity_decoy_start)
        
        ####
        spe_bf_num_decoy_start<-sum(spe_ref_XCorr_bf_decoy_start!=0)
        
        if(max(spe_ref_XCorr_bf_decoy_start)>0){
          max_cor_index_start<-which.max(spe_ref_XCorr_bf_decoy_start)
          if(max_cor_index_start<=(nfi_bfmass_start-nfimass)/3){
            max_cor_by_ph_bf_decoy_start<-blabel1[nphk+max_cor_index_start-1]
          }
          else if(max_cor_index_start<=2*(nfi_bfmass_start-nfimass)/3 & max_cor_index_start>(nfi_bfmass_start-nfimass)/3){
            max_cor_index1_start<-max_cor_index_start-(nfi_bfmass_start-nfimass)/3
            max_cor_by_ph_bf_decoy_start<-blabel2[nphk+max_cor_index1_start-1]
          }
          else{
            max_cor_index2_start<-max_cor_index_start-2*(nfi_bfmass_start-nfimass)/3
            max_cor_by_ph_bf_decoy_start<-blabel3[nphk+max_cor_index2_start-1]
          }
        }else{
          max_cor_by_ph_bf_decoy_start<-"NA"
        }
        
        ####
        ###
        fimass_ph1_intensity_start<-vector()
        str_fi_intensity_start<-vector()
        fimass_ph1_noise_start<-vector()
        str_fi_noise_start<-vector()
        ####
        fimass_ph1_intensity_decoy_start<-vector()
        str_fi_intensity_decoy_start<-vector()
        fimass_ph1_noise_decoy_start<-vector()
        str_fi_noise_decoy_start<-vector()
        
        for(sfi1_start in 1:nfi_yhmass_start){
          for(nrt_start in 1:nnewrt){
            nrt1_start<-newrt[nrt_start]
            char_mz_start<-as.character(datapsm01$m.z[nrt1_start])
            str_mz_start<-strsplit(char_mz_start,";")
            newstr_mz_start<-as.numeric(str_mz_start[[1]])
            char_intensity_start<-as.character(datapsm01$Intensity[nrt1_start])
            str_intensity_start<-strsplit(char_intensity_start,";")
            newstr_intensity_start<-as.numeric(str_intensity_start[[1]])
            char_noise_start<-as.character(datapsm01$Noise[nrt1_start])
            str_noise_start<-strsplit(char_noise_start,";")
            newstr_noise_start<-as.numeric(str_noise_start[[1]])
            
            sub_fion1_start<-abs(newstr_mz_start-fi_yhmass_start[sfi1_start])
            if(min(sub_fion1_start)<=peakTolThreshold){
              sub_fion1_min_start<-which.min(sub_fion1_start)
              fimass_ph1_intensity_start[nrt_start]<-newstr_intensity_start[sub_fion1_min_start]
              fimass_ph1_noise_start[nrt_start]<-newstr_noise_start[sub_fion1_min_start]
            }else{
              fimass_ph1_intensity_start[nrt_start]<-0
              fimass_ph1_noise_start[nrt_start]<-0.0001
            }
            ####
            sub_fion1_decoy_start<-abs(newstr_mz_start-fi_yhmass_decoy_start[sfi1_start])
            if(min(sub_fion1_decoy_start)<=peakTolThreshold){
              sub_fion1_min_decoy_start<-which.min(sub_fion1_decoy_start)
              fimass_ph1_intensity_decoy_start[nrt_start]<-newstr_intensity_start[sub_fion1_min_decoy_start]
              fimass_ph1_noise_decoy_start[nrt_start]<-newstr_noise_start[sub_fion1_min_decoy_start]
            }else{
              fimass_ph1_intensity_decoy_start[nrt_start]<-0
              fimass_ph1_noise_decoy_start[nrt_start]<-0.0001
            }
          }
          newfimass_ph1_inten_start<-paste(fimass_ph1_intensity_start,collapse=";")
          str_fi_intensity_start[sfi1_start]<-newfimass_ph1_inten_start
          newfimass_ph1_noise_start<-paste(fimass_ph1_noise_start,collapse=";")
          str_fi_noise_start[sfi1_start]<-newfimass_ph1_noise_start
          ####
          newfimass_ph1_inten_decoy_start<-paste(fimass_ph1_intensity_decoy_start,collapse=";")
          str_fi_intensity_decoy_start[sfi1_start]<-newfimass_ph1_inten_decoy_start
          newfimass_ph1_noise_decoy_start<-paste(fimass_ph1_noise_decoy_start,collapse=";")
          str_fi_noise_decoy_start[sfi1_start]<-newfimass_ph1_noise_decoy_start
        }
        
        #####Target_y ion ahead
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yh_start<-vector()
        spe_yh_intensity_start<-vector()
        ncor_yh_start<-1
        for(npyh2_start in (nfimass+1):nfi_yhmass_start){
          str_fi_intensitynew2_start<-strsplit(str_fi_intensity_start[npyh2_start],";")
          str_fi_intensitynew3_start<-as.numeric(str_fi_intensitynew2_start[[1]])
          nnozero_str_fi_intensitynew3_1_start<-sum(str_fi_intensitynew3_start!=0)
          str_fi_noisenew2_start<-strsplit(str_fi_noise_start[npyh2_start],";")
          str_fi_noisenew3_start<-as.numeric(str_fi_noisenew2_start[[1]])
          if(nnozero_str_fi_intensitynew3_1_start<ppt){
            spe_ref_XCorr_yh_start[ncor_yh_start]<-0
            spe_yh_intensity_start[ncor_yh_start]<-0
          }else{
            cor_yh1_start<-vector()
            sum_str_fi_intensitynew1_start<-vector()
            for(npyh1_start in 1:nfimass){
              str_fi_intensitynew_start<-strsplit(str_fi_intensity_start[npyh1_start],";")
              str_fi_intensitynew1_start<-as.numeric(str_fi_intensitynew_start[[1]])
              nnozero_str_fi_intensitynew1_start<-sum(str_fi_intensitynew1_start!=0)
              if(nnozero_str_fi_intensitynew1_start<ppt){
                cor_yh1_start[npyh1_start]<-0
                sum_str_fi_intensitynew1_start[npyh1_start]<-0
              }else{
                cor_yh1_start[npyh1_start]<-sum(str_fi_intensitynew1_start*str_fi_intensitynew3_start)/sqrt(sum(str_fi_intensitynew1_start^2)*sum(str_fi_intensitynew3_start^2))
                sum_str_fi_intensitynew1_start[npyh1_start]<-log2(sum(str_fi_intensitynew1_start))
              }
            }
            if(all(cor_yh1_start==0)){
              spe_ref_XCorr_yh_start[ncor_yh_start]<-0
              spe_yh_intensity_start[ncor_yh_start]<-0
            }else{
              cor_yh_start_weight<-sum(cor_yh1_start)#*sum_str_fi_intensitynew1_start)/sum(sum_str_fi_intensitynew1_start)
              #cor_yh_noise_start<-str_fi_intensitynew3_start/str_fi_noisenew3_start
              #ncor_yh_noise_start<-sum(cor_yh_noise_start>0)
              #cor_yh_noise_index_start<-which(cor_yh_noise_start>0)
              #cor_yh_noise_sum_start<-sum(cor_yh_noise_start[cor_yh_noise_index_start])
              if(cor_yh_start_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yh_start[ncor_yh_start]<-cor_yh_start_weight
                nozero_str_fi_intensitynew3_start<-log2(sum(str_fi_intensitynew3_start))#*nnozero_str_fi_intensitynew3_start
                spe_yh_intensity_start[ncor_yh_start]<-nozero_str_fi_intensitynew3_start
              }else{
                spe_ref_XCorr_yh_start[ncor_yh_start]<-0
                spe_yh_intensity_start[ncor_yh_start]<-0
              }
              #*(1-exp(-(cor_yh_noise_sum_start-ncor_yh_noise_start)))
            }
          }
          ncor_yh_start<-ncor_yh_start+1
        }
        spe_ref_XCorr_yh_sum_start<-sum(spe_ref_XCorr_yh_start)
        spe_yh_intensity_sum_start<-sum(spe_yh_intensity_start)
        
        ####
        spe_yh_num_start<-sum(spe_ref_XCorr_yh_start!=0)
        
        if(max(spe_ref_XCorr_yh_start)>0){
          max_y_cor_index_start<-which.max(spe_ref_XCorr_yh_start)
          if(max_y_cor_index_start<=(nfi_yhmass_start-nfimass)/3){
            newmax_cor_by_start<-ylabel1[nby-nphk-max_y_cor_index_start+1]
          }
          else if(max_y_cor_index_start<=2*(nfi_yhmass_start-nfimass)/3 & max_y_cor_index_start>(nfi_yhmass_start-nfimass)/3){
            max_y_cor_index1_start<-max_y_cor_index_start-(nfi_yhmass_start-nfimass)/3
            newmax_cor_by_start<-ylabel2[nby-nphk-max_y_cor_index1_start+1]
          }
          else{
            max_y_cor_index2_start<-max_y_cor_index_start-2*(nfi_yhmass_start-nfimass)/3
            newmax_cor_by_start<-ylabel3[nby-nphk-max_y_cor_index2_start+1]
          }
        }else{
          newmax_cor_by_start<-"NA"
        }
        
        #####Decoy_y ion ahead
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yh_decoy_start<-vector()
        spe_yh_intensity_decoy_start<-vector()
        ncor_yh_decoy_start<-1
        for(npyh2_start in (nfimass+1):nfi_yhmass_start){
          str_fi_intensitynew2_start<-strsplit(str_fi_intensity_decoy_start[npyh2_start],";")
          str_fi_intensitynew3_start<-as.numeric(str_fi_intensitynew2_start[[1]])
          nnozero_str_fi_intensitynew3_1_start<-sum(str_fi_intensitynew3_start!=0)
          str_fi_noisenew2_start<-strsplit(str_fi_noise_decoy_start[npyh2_start],";")
          str_fi_noisenew3_start<-as.numeric(str_fi_noisenew2_start[[1]])
          if(nnozero_str_fi_intensitynew3_1_start<ppt){
            spe_ref_XCorr_yh_decoy_start[ncor_yh_decoy_start]<-0
            spe_yh_intensity_decoy_start[ncor_yh_decoy_start]<-0
          }else{
            cor_yh1_start<-vector()
            sum_str_fi_intensitynew1_start<-vector()
            for(npyh1_start in 1:nfimass){
              str_fi_intensitynew_start<-strsplit(str_fi_intensity_decoy_start[npyh1_start],";")
              str_fi_intensitynew1_start<-as.numeric(str_fi_intensitynew_start[[1]])
              nnozero_str_fi_intensitynew1_start<-sum(str_fi_intensitynew1_start!=0)
              if(nnozero_str_fi_intensitynew1_start<ppt){
                cor_yh1_start[npyh1_start]<-0
                sum_str_fi_intensitynew1_start[npyh1_start]<-0
              }else{
                cor_yh1_start[npyh1_start]<-sum(str_fi_intensitynew1_start*str_fi_intensitynew3_start)/sqrt(sum(str_fi_intensitynew1_start^2)*sum(str_fi_intensitynew3_start^2))
                sum_str_fi_intensitynew1_start[npyh1_start]<-log2(sum(str_fi_intensitynew1_start))
              }
            }
            if(all(cor_yh1_start==0)){
              spe_ref_XCorr_yh_decoy_start[ncor_yh_decoy_start]<-0
              spe_yh_intensity_decoy_start[ncor_yh_decoy_start]<-0
            }else{
              cor_yh_start_weight<-sum(cor_yh1_start)#*sum_str_fi_intensitynew1_start)/sum(sum_str_fi_intensitynew1_start)
              #cor_yh_noise_start<-str_fi_intensitynew3_start/str_fi_noisenew3_start
              #ncor_yh_noise_start<-sum(cor_yh_noise_start>0)
              #cor_yh_noise_index_start<-which(cor_yh_noise_start>0)
              #cor_yh_noise_sum_start<-sum(cor_yh_noise_start[cor_yh_noise_index_start])
              if(cor_yh_start_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yh_decoy_start[ncor_yh_decoy_start]<-cor_yh_start_weight
                nozero_str_fi_intensitynew3_start<-log2(sum(str_fi_intensitynew3_start))#*nnozero_str_fi_intensitynew3_start
                spe_yh_intensity_decoy_start[ncor_yh_decoy_start]<-nozero_str_fi_intensitynew3_start
              }else{
                spe_ref_XCorr_yh_decoy_start[ncor_yh_decoy_start]<-0
                spe_yh_intensity_decoy_start[ncor_yh_decoy_start]<-0
              }
              #*(1-exp(-(cor_yh_noise_sum_start-ncor_yh_noise_start)))
              
            }
          }
          ncor_yh_decoy_start<-ncor_yh_decoy_start+1
        }
        spe_ref_XCorr_yh_sum_decoy_start<-sum(spe_ref_XCorr_yh_decoy_start)
        spe_yh_intensity_sum_decoy_start<-sum(spe_yh_intensity_decoy_start)
        
        ####
        spe_yh_num_decoy_start<-sum(spe_ref_XCorr_yh_decoy_start!=0)
        
        if(max(spe_ref_XCorr_yh_decoy_start)>0){
          max_y_cor_index_start<-which.max(spe_ref_XCorr_yh_decoy_start)
          if(max_y_cor_index_start<=(nfi_yhmass_start-nfimass)/3){
            newmax_cor_by_decoy_start<-ylabel1[nby-nphk-max_y_cor_index_start+1]
          }
          else if(max_y_cor_index_start<=2*(nfi_yhmass_start-nfimass)/3 & max_y_cor_index_start>(nfi_yhmass_start-nfimass)/3){
            max_y_cor_index1_start<-max_y_cor_index_start-(nfi_yhmass_start-nfimass)/3
            newmax_cor_by_decoy_start<-ylabel2[nby-nphk-max_y_cor_index1_start+1]
          }
          else{
            max_y_cor_index2_start<-max_y_cor_index_start-2*(nfi_yhmass_start-nfimass)/3
            newmax_cor_by_decoy_start<-ylabel3[nby-nphk-max_y_cor_index2_start+1]
          }
        }else{
          newmax_cor_by_decoy_start<-"NA"
        }
        
        ####
        spe_ref_XCorr_sum_ph[phk]<-spe_ref_XCorr_bf_sum_start+spe_ref_XCorr_yh_sum_start
        spe_fi_num_ph[phk]<-spe_bf_num_start+spe_yh_num_start
        spe_fi_intensity_sum_ph[phk]<-spe_bf_intensity_sum_start+spe_yh_intensity_sum_start
        max_cor_by_ph[phk]<-paste(max_cor_by_ph_bf_start,newmax_cor_by_start,sep=" ")
        ####
        spe_ref_XCorr_sum_ph_decoy[phk]<-spe_ref_XCorr_bf_sum_decoy_start+spe_ref_XCorr_yh_sum_decoy_start
        spe_fi_num_ph_decoy[phk]<-spe_bf_num_decoy_start+spe_yh_num_decoy_start
        spe_fi_intensity_sum_ph_decoy[phk]<-spe_bf_intensity_sum_decoy_start+spe_yh_intensity_sum_decoy_start
        max_cor_by_ph_decoy[phk]<-paste(max_cor_by_ph_bf_decoy_start,newmax_cor_by_decoy_start,sep=" ")
      }
      else if(nphk==char_sty[[1]][nchar_sty]){
        phk_bmass_ahead11_end<-newbmass[char_sty[[1]][nchar_sty-1]:(nphk-1)]
        phk_bmass_ahead12_end<-phk_bmass_ahead11_end-18.0152
        phk_bmass_ahead13_end<-phk_bmass_ahead11_end-17.0304
        phk_bmass_ahead1_end<-c(phk_bmass_ahead11_end,phk_bmass_ahead12_end,
                                phk_bmass_ahead13_end)
        phk_ymass_after11_end<-yyrev[(char_sty[[1]][nchar_sty-1]+1):nphk]
        phk_ymass_after12_end<-phk_ymass_after11_end-18.0152
        phk_ymass_after13_end<-phk_ymass_after11_end-17.0304
        phk_ymass_after1_end<-c(phk_ymass_after11_end,phk_ymass_after12_end,
                                phk_ymass_after13_end)
        fi_bhmass_end<-c(fimass,phk_bmass_ahead1_end)
        fi_yfmass_end<-c(fimass,phk_ymass_after1_end)
        nfi_bhmass_end<-length(fi_bhmass_end)
        nfi_yfmass_end<-length(fi_yfmass_end)
        
        fimass_ph1_intensity_end<-vector()
        str_fi_intensity_end<-vector()
        fimass_ph1_noise_end<-vector()
        str_fi_noise_end<-vector()
        ####
        phk_bmass_ahead_decoy11_end<-newbmass_decoy[char_sty[[1]][nchar_sty-1]:(nphk-1)]
        phk_bmass_ahead_decoy12_end<-phk_bmass_ahead_decoy11_end-18.0152
        phk_bmass_ahead_decoy13_end<-phk_bmass_ahead_decoy11_end-17.0304
        phk_bmass_ahead_decoy1_end<-c(phk_bmass_ahead_decoy11_end,phk_bmass_ahead_decoy12_end,
                                      phk_bmass_ahead_decoy13_end)
        phk_ymass_after_decoy11_end<-yyrev_decoy[(char_sty[[1]][nchar_sty-1]+1):nphk]
        phk_ymass_after_decoy12_end<-phk_ymass_after_decoy11_end-18.0152
        phk_ymass_after_decoy13_end<-phk_ymass_after_decoy11_end-17.0304
        phk_ymass_after_decoy1_end<-c(phk_ymass_after_decoy11_end,phk_ymass_after_decoy12_end,
                                      phk_ymass_after_decoy13_end)
        fi_bhmass_decoy_end<-c(fimass_decoy,phk_bmass_ahead_decoy1_end)
        fi_yfmass_decoy_end<-c(fimass_decoy,phk_ymass_after_decoy1_end)
        nfi_bhmass_decoy_end<-length(fi_bhmass_decoy_end)
        nfi_yfmass_decoy_end<-length(fi_yfmass_decoy_end)
        
        fimass_ph1_intensity_decoy_end<-vector()
        str_fi_intensity_decoy_end<-vector()
        fimass_ph1_noise_decoy_end<-vector()
        str_fi_noise_decoy_end<-vector()
        ####
        for(sfi1_end in 1:nfi_bhmass_end){
          for(nrt_end in 1:nnewrt){
            nrt1_end<-newrt[nrt_end]
            char_mz_end<-as.character(datapsm01$m.z[nrt1_end])
            str_mz_end<-strsplit(char_mz_end,";")
            newstr_mz_end<-as.numeric(str_mz_end[[1]])
            char_intensity_end<-as.character(datapsm01$Intensity[nrt1_end])
            str_intensity_end<-strsplit(char_intensity_end,";")
            newstr_intensity_end<-as.numeric(str_intensity_end[[1]])
            char_noise_end<-as.character(datapsm01$Noise[nrt1_end])
            str_noise_end<-strsplit(char_noise_end,";")
            newstr_noise_end<-as.numeric(str_noise_end[[1]])
            
            sub_fion1_end<-abs(newstr_mz_end-fi_bhmass_end[sfi1_end])
            if(min(sub_fion1_end)<=peakTolThreshold){
              sub_fion1_min_end<-which.min(sub_fion1_end)
              fimass_ph1_intensity_end[nrt_end]<-newstr_intensity_end[sub_fion1_min_end]
              fimass_ph1_noise_end[nrt_end]<-newstr_noise_end[sub_fion1_min_end]
            }else{
              fimass_ph1_intensity_end[nrt_end]<-0
              fimass_ph1_noise_end[nrt_end]<-0.0001
            }
            ####
            sub_fion1_decoy_end<-abs(newstr_mz_end-fi_bhmass_decoy_end[sfi1_end])
            if(min(sub_fion1_decoy_end)<=peakTolThreshold){
              sub_fion1_min_decoy_end<-which.min(sub_fion1_decoy_end)
              fimass_ph1_intensity_decoy_end[nrt_end]<-newstr_intensity_end[sub_fion1_min_decoy_end]
              fimass_ph1_noise_decoy_end[nrt_end]<-newstr_noise_end[sub_fion1_min_decoy_end]
            }else{
              fimass_ph1_intensity_decoy_end[nrt_end]<-0
              fimass_ph1_noise_decoy_end[nrt_end]<-0.0001
            }
          }
          newfimass_ph1_inten_end<-paste(fimass_ph1_intensity_end,collapse=";")
          str_fi_intensity_end[sfi1_end]<-newfimass_ph1_inten_end
          newfimass_ph1_noise_end<-paste(fimass_ph1_noise_end,collapse=";")
          str_fi_noise_end[sfi1_end]<-newfimass_ph1_noise_end
          ####
          newfimass_ph1_inten_decoy_end<-paste(fimass_ph1_intensity_decoy_end,collapse=";")
          str_fi_intensity_decoy_end[sfi1_end]<-newfimass_ph1_inten_decoy_end
          newfimass_ph1_noise_decoy_end<-paste(fimass_ph1_noise_decoy_end,collapse=";")
          str_fi_noise_decoy_end[sfi1_end]<-newfimass_ph1_noise_decoy_end
        }
        
        #####Target_b ion ahead
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bh_end<-vector()
        spe_bh_intensity_end<-vector()
        ncor_bh_end<-1
        for(npbh2_end in (nfimass+1):nfi_bhmass_end){
          str_fi_intensitynew2_end<-strsplit(str_fi_intensity_end[npbh2_end],";")
          str_fi_intensitynew3_end<-as.numeric(str_fi_intensitynew2_end[[1]])
          nnozero_str_fi_intensitynew3_1_end<-sum(str_fi_intensitynew3_end!=0)
          str_fi_noisenew2_end<-strsplit(str_fi_noise_end[npbh2_end],";")
          str_fi_noisenew3_end<-as.numeric(str_fi_noisenew2_end[[1]])
          if(nnozero_str_fi_intensitynew3_1_end<ppt){
            spe_ref_XCorr_bh_end[ncor_bh_end]<-0
            spe_bh_intensity_end[ncor_bh_end]<-0
          }else{
            cor_bh1_end<-vector()
            sum_str_fi_intensitynew1_end<-vector()
            for(npbh1_end in 1:nfimass){
              str_fi_intensitynew_end<-strsplit(str_fi_intensity_end[npbh1_end],";")
              str_fi_intensitynew1_end<-as.numeric(str_fi_intensitynew_end[[1]])
              nnozero_str_fi_intensitynew1_end<-sum(str_fi_intensitynew1_end!=0)
              if(nnozero_str_fi_intensitynew1_end<ppt){
                cor_bh1_end[npbh1_end]<-0
                sum_str_fi_intensitynew1_end[npbh1_end]<-0
              }else{
                cor_bh1_end[npbh1_end]<-sum(str_fi_intensitynew1_end*str_fi_intensitynew3_end)/sqrt(sum(str_fi_intensitynew1_end^2)*sum(str_fi_intensitynew3_end^2))
                sum_str_fi_intensitynew1_end[npbh1_end]<-log2(sum(str_fi_intensitynew1_end))
              }
            }
            if(all(cor_bh1_end==0)){
              spe_ref_XCorr_bh_end[ncor_bh_end]<-0
              spe_bh_intensity_end[ncor_bh_end]<-0
            }else{
              cor_bh_end_weight<-sum(cor_bh1_end)#*sum_str_fi_intensitynew1_end)/sum(sum_str_fi_intensitynew1_end)
              #cor_bh_noise_end<-str_fi_intensitynew3_end/str_fi_noisenew3_end
              #ncor_bh_noise_end<-sum(cor_bh_noise_end>0)
              #cor_bh_noise_index_end<-which(cor_bh_noise_end>0)
              #cor_bh_noise_sum_end<-sum(cor_bh_noise_end[cor_bh_noise_index_end])
              if(cor_bh_end_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bh_end[ncor_bh_end]<-cor_bh_end_weight
                nozero_str_fi_intensitynew3_end<-log2(sum(str_fi_intensitynew3_end))#*nnozero_str_fi_intensitynew3_end
                spe_bh_intensity_end[ncor_bh_end]<-nozero_str_fi_intensitynew3_end
              }else{
                spe_ref_XCorr_bh_end[ncor_bh_end]<-0
                spe_bh_intensity_end[ncor_bh_end]<-0
              }
              #*(1-exp(-(cor_bh_noise_sum_end-ncor_bh_noise_end)))
              
            }
          }
          ncor_bh_end<-ncor_bh_end+1
        }
        spe_ref_XCorr_bh_sum_end<-sum(spe_ref_XCorr_bh_end)
        spe_bh_intensity_sum_end<-sum(spe_bh_intensity_end)
        
        ####
        spe_bh_num_end<-sum(spe_ref_XCorr_bh_end!=0)
        
        if(max(spe_ref_XCorr_bh_end)>0){
          max_cor_index_end<-which.max(spe_ref_XCorr_bh_end)
          if(max_cor_index_end<=(nfi_bhmass_end-nfimass)/3){
            max_cor_by_ph_bh_end<-blabel1[nphk+max_cor_index_end-(nfi_bhmass_end-nfimass)/3-1]
          }
          else if(max_cor_index_end<=2*(nfi_bhmass_end-nfimass)/3 & max_cor_index_end>(nfi_bhmass_end-nfimass)/3){
            max_cor_index1_end<-max_cor_index_end-(nfi_bhmass_end-nfimass)/3
            max_cor_by_ph_bh_end<-blabel2[nphk+max_cor_index1_end-(nfi_bhmass_end-nfimass)/3-1]
          }
          else{
            max_cor_index2_end<-max_cor_index_end-2*(nfi_bhmass_end-nfimass)/3
            max_cor_by_ph_bh_end<-blabel3[nphk+max_cor_index2_end-(nfi_bhmass_end-nfimass)/3-1]
          }
        }else{
          max_cor_by_ph_bh_end<-"NA"
        }
        
        #####Decoy_b ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bh_decoy_end<-vector()
        spe_bh_intensity_decoy_end<-vector()
        ncor_bh_decoy_end<-1
        for(npbh2_end in (nfimass+1):nfi_bhmass_end){
          str_fi_intensitynew2_end<-strsplit(str_fi_intensity_decoy_end[npbh2_end],";")
          str_fi_intensitynew3_end<-as.numeric(str_fi_intensitynew2_end[[1]])
          nnozero_str_fi_intensitynew3_1_end<-sum(str_fi_intensitynew3_end!=0)
          str_fi_noisenew2_end<-strsplit(str_fi_noise_decoy_end[npbh2_end],";")
          str_fi_noisenew3_end<-as.numeric(str_fi_noisenew2_end[[1]])
          if(nnozero_str_fi_intensitynew3_1_end<ppt){
            spe_ref_XCorr_bh_decoy_end[ncor_bh_decoy_end]<-0
            spe_bh_intensity_decoy_end[ncor_bh_decoy_end]<-0
          }else{
            cor_bh1_end<-vector()
            sum_str_fi_intensitynew1_end<-vector()
            for(npbh1_end in 1:nfimass){
              str_fi_intensitynew_end<-strsplit(str_fi_intensity_decoy_end[npbh1_end],";")
              str_fi_intensitynew1_end<-as.numeric(str_fi_intensitynew_end[[1]])
              nnozero_str_fi_intensitynew1_end<-sum(str_fi_intensitynew1_end!=0)
              if(nnozero_str_fi_intensitynew1_end<ppt){
                cor_bh1_end[npbh1_end]<-0
                sum_str_fi_intensitynew1_end[npbh1_end]<-0
              }else{
                cor_bh1_end[npbh1_end]<-sum(str_fi_intensitynew1_end*str_fi_intensitynew3_end)/sqrt(sum(str_fi_intensitynew1_end^2)*sum(str_fi_intensitynew3_end^2))
                sum_str_fi_intensitynew1_end[npbh1_end]<-log2(sum(str_fi_intensitynew1_end))
              }
            }
            if(all(cor_bh1_end==0)){
              spe_ref_XCorr_bh_decoy_end[ncor_bh_decoy_end]<-0
              spe_bh_intensity_decoy_end[ncor_bh_decoy_end]<-0
            }else{
              cor_bh_end_weight<-sum(cor_bh1_end)#*sum_str_fi_intensitynew1_end)/sum(sum_str_fi_intensitynew1_end)
              #cor_bh_noise_end<-str_fi_intensitynew3_end/str_fi_noisenew3_end
              #ncor_bh_noise_end<-sum(cor_bh_noise_end>0)
              #cor_bh_noise_index_end<-which(cor_bh_noise_end>0)
              #cor_bh_noise_sum_end<-sum(cor_bh_noise_end[cor_bh_noise_index_end])
              if(cor_bh_end_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bh_decoy_end[ncor_bh_decoy_end]<-cor_bh_end_weight
                nozero_str_fi_intensitynew3_end<-log2(sum(str_fi_intensitynew3_end))#*nnozero_str_fi_intensitynew3_end
                spe_bh_intensity_decoy_end[ncor_bh_decoy_end]<-nozero_str_fi_intensitynew3_end
              }else{
                spe_ref_XCorr_bh_decoy_end[ncor_bh_decoy_end]<-0
                spe_bh_intensity_decoy_end[ncor_bh_decoy_end]<-0
              }
              #*(1-exp(-(cor_bh_noise_sum_end-ncor_bh_noise_end)))
              
            }
          }
          ncor_bh_decoy_end<-ncor_bh_decoy_end+1
        }
        spe_ref_XCorr_bh_sum_decoy_end<-sum(spe_ref_XCorr_bh_decoy_end)
        spe_bh_intensity_sum_decoy_end<-sum(spe_bh_intensity_decoy_end)
        
        ####
        spe_bh_num_decoy_end<-sum(spe_ref_XCorr_bh_decoy_end!=0)
        
        if(max(spe_ref_XCorr_bh_decoy_end)>0){
          max_cor_index_end<-which.max(spe_ref_XCorr_bh_decoy_end)
          if(max_cor_index_end<=(nfi_bhmass_end-nfimass)/3){
            max_cor_by_ph_bh_decoy_end<-blabel1[nphk+max_cor_index_end-(nfi_bhmass_end-nfimass)/3-1]
          }
          else if(max_cor_index_end<=2*(nfi_bhmass_end-nfimass)/3 & max_cor_index_end>(nfi_bhmass_end-nfimass)/3){
            max_cor_index1_end<-max_cor_index_end-(nfi_bhmass_end-nfimass)/3
            max_cor_by_ph_bh_decoy_end<-blabel2[nphk+max_cor_index1_end-(nfi_bhmass_end-nfimass)/3-1]
          }
          else{
            max_cor_index3_end<-max_cor_index_end-2*(nfi_bhmass_end-nfimass)/3
            max_cor_by_ph_bh_decoy_end<-blabel3[nphk+max_cor_index3_end-(nfi_bhmass_end-nfimass)/3-1]
          }
        }else{
          max_cor_by_ph_bh_decoy_end<-"NA"
        }
        
        ####
        ###
        fimass_ph1_intensity_end<-vector()
        str_fi_intensity_end<-vector()
        fimass_ph1_noise_end<-vector()
        str_fi_noise_end<-vector()
        ####
        fimass_ph1_intensity_decoy_end<-vector()
        str_fi_intensity_decoy_end<-vector()
        fimass_ph1_noise_decoy_end<-vector()
        str_fi_noise_decoy_end<-vector()
        
        for(sfi1_end in 1:nfi_yfmass_end){
          for(nrt_end in 1:nnewrt){
            nrt1_end<-newrt[nrt_end]
            char_mz_end<-as.character(datapsm01$m.z[nrt1_end])
            str_mz_end<-strsplit(char_mz_end,";")
            newstr_mz_end<-as.numeric(str_mz_end[[1]])
            char_intensity_end<-as.character(datapsm01$Intensity[nrt1_end])
            str_intensity_end<-strsplit(char_intensity_end,";")
            newstr_intensity_end<-as.numeric(str_intensity_end[[1]])
            char_noise_end<-as.character(datapsm01$Noise[nrt1_end])
            str_noise_end<-strsplit(char_noise_end,";")
            newstr_noise_end<-as.numeric(str_noise_end[[1]])
            
            sub_fion1_end<-abs(newstr_mz_end-fi_yfmass_end[sfi1_end])
            if(min(sub_fion1_end)<=peakTolThreshold){
              sub_fion1_min_end<-which.min(sub_fion1_end)
              fimass_ph1_intensity_end[nrt_end]<-newstr_intensity_end[sub_fion1_min_end]
              fimass_ph1_noise_end[nrt_end]<-newstr_noise_end[sub_fion1_min_end]
            }else{
              fimass_ph1_intensity_end[nrt_end]<-0
              fimass_ph1_noise_end[nrt_end]<-0.0001
            }
            ####
            sub_fion1_decoy_end<-abs(newstr_mz_end-fi_yfmass_decoy_end[sfi1_end])
            if(min(sub_fion1_decoy_end)<=peakTolThreshold){
              sub_fion1_min_decoy_end<-which.min(sub_fion1_decoy_end)
              fimass_ph1_intensity_decoy_end[nrt_end]<-newstr_intensity_end[sub_fion1_min_decoy_end]
              fimass_ph1_noise_decoy_end[nrt_end]<-newstr_noise_end[sub_fion1_min_decoy_end]
            }else{
              fimass_ph1_intensity_decoy_end[nrt_end]<-0
              fimass_ph1_noise_decoy_end[nrt_end]<-0.0001
            }
          }
          newfimass_ph1_inten_end<-paste(fimass_ph1_intensity_end,collapse=";")
          str_fi_intensity_end[sfi1_end]<-newfimass_ph1_inten_end
          newfimass_ph1_noise_end<-paste(fimass_ph1_noise_end,collapse=";")
          str_fi_noise_end[sfi1_end]<-newfimass_ph1_noise_end
          ####
          newfimass_ph1_inten_decoy_end<-paste(fimass_ph1_intensity_decoy_end,collapse=";")
          str_fi_intensity_decoy_end[sfi1_end]<-newfimass_ph1_inten_decoy_end
          newfimass_ph1_noise_decoy_end<-paste(fimass_ph1_noise_decoy_end,collapse=";")
          str_fi_noise_decoy_end[sfi1_end]<-newfimass_ph1_noise_decoy_end
        }
        
        #####Target_y ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yf_end<-vector()
        spe_yf_intensity_end<-vector()
        ncor_yf_end<-1
        for(npyf2_end in (nfimass+1):nfi_yfmass_end){
          str_fi_intensitynew2_end<-strsplit(str_fi_intensity_end[npyf2_end],";")
          str_fi_intensitynew3_end<-as.numeric(str_fi_intensitynew2_end[[1]])
          nnozero_str_fi_intensitynew3_1_end<-sum(str_fi_intensitynew3_end!=0)
          str_fi_noisenew2_end<-strsplit(str_fi_noise_end[npyf2_end],";")
          str_fi_noisenew3_end<-as.numeric(str_fi_noisenew2_end[[1]])
          if(nnozero_str_fi_intensitynew3_1_end<ppt){
            spe_ref_XCorr_yf_end[ncor_yf_end]<-0
            spe_yf_intensity_end[ncor_yf_end]<-0
          }else{
            cor_yf1_end<-vector()
            sum_str_fi_intensitynew1_end<-vector()
            for(npyf1_end in 1:nfimass){
              str_fi_intensitynew_end<-strsplit(str_fi_intensity_end[npyf1_end],";")
              str_fi_intensitynew1_end<-as.numeric(str_fi_intensitynew_end[[1]])
              nnozero_str_fi_intensitynew1_end<-sum(str_fi_intensitynew1_end!=0)
              if(nnozero_str_fi_intensitynew1_end<ppt){
                cor_yf1_end[npyf1_end]<-0
                sum_str_fi_intensitynew1_end[npyf1_end]<-0
              }else{
                cor_yf1_end[npyf1_end]<-sum(str_fi_intensitynew1_end*str_fi_intensitynew3_end)/sqrt(sum(str_fi_intensitynew1_end^2)*sum(str_fi_intensitynew3_end^2))
                sum_str_fi_intensitynew1_end[npyf1_end]<-log2(sum(str_fi_intensitynew1_end))
              }
            }
            if(all(cor_yf1_end==0)){
              spe_ref_XCorr_yf_end[ncor_yf_end]<-0
              spe_yf_intensity_end[ncor_yf_end]<-0
            }else{
              cor_yf_end_weight<-sum(cor_yf1_end)#*sum_str_fi_intensitynew1_end)/sum(sum_str_fi_intensitynew1_end)
              #cor_yf_noise_end<-str_fi_intensitynew3_end/str_fi_noisenew3_end
              #ncor_yf_noise_end<-sum(cor_yf_noise_end>0)
              #cor_yf_noise_index_end<-which(cor_yf_noise_end>0)
              #cor_yf_noise_sum_end<-sum(cor_yf_noise_end[cor_yf_noise_index_end])
              if(cor_yf_end_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yf_end[ncor_yf_end]<-cor_yf_end_weight
                nozero_str_fi_intensitynew3_end<-log2(sum(str_fi_intensitynew3_end))#*nnozero_str_fi_intensitynew3_end
                spe_yf_intensity_end[ncor_yf_end]<-nozero_str_fi_intensitynew3_end
              }else{
                spe_ref_XCorr_yf_end[ncor_yf_end]<-0
                spe_yf_intensity_end[ncor_yf_end]<-0
              }
              #*(1-exp(-(cor_yf_noise_sum_end-ncor_yf_noise_end)))
              
            }
          }
          ncor_yf_end<-ncor_yf_end+1
        }
        spe_ref_XCorr_yf_sum_end<-sum(spe_ref_XCorr_yf_end)
        spe_yf_intensity_sum_end<-sum(spe_yf_intensity_end)
        
        ####
        spe_yf_num_end<-sum(spe_ref_XCorr_yf_end!=0)
        
        if(max(spe_ref_XCorr_yf_end)>0){
          max_y_cor_index_end<-which.max(spe_ref_XCorr_yf_end)
          if(max_y_cor_index_end<=(nfi_yfmass_end-nfimass)/3){
            newmax_cor_by_end<-ylabel1[nby-nphk+(nfi_yfmass_end-nfimass)/3-max_y_cor_index_end+1]
          }
          else if(max_y_cor_index_end<=2*(nfi_yfmass_end-nfimass)/3 & max_y_cor_index_end>(nfi_yfmass_end-nfimass)/3){
            max_y_cor_index1_end<-max_y_cor_index_end-(nfi_yfmass_end-nfimass)/3
            newmax_cor_by_end<-ylabel2[nby-nphk+(nfi_yfmass_end-nfimass)/3-max_y_cor_index1_end+1]
          }
          else{
            max_y_cor_index3_end<-max_y_cor_index_end-2*nfi_yfmass_end-nfimass/3
            newmax_cor_by_end<-ylabel3[nby-nphk+(nfi_yfmass_end-nfimass)/3-max_y_cor_index3_end+1]
          }
        }else{
          newmax_cor_by_end<-"NA"
        }
        
        #####Decoy_y ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yf_decoy_end<-vector()
        spe_yf_intensity_decoy_end<-vector()
        ncor_yf_decoy_end<-1
        for(npyf2_end in (nfimass+1):nfi_yfmass_end){
          str_fi_intensitynew2_end<-strsplit(str_fi_intensity_decoy_end[npyf2_end],";")
          str_fi_intensitynew3_end<-as.numeric(str_fi_intensitynew2_end[[1]])
          nnozero_str_fi_intensitynew3_1_end<-sum(str_fi_intensitynew3_end!=0)
          str_fi_noisenew2_end<-strsplit(str_fi_noise_decoy_end[npyf2_end],";")
          str_fi_noisenew3_end<-as.numeric(str_fi_noisenew2_end[[1]])
          if(nnozero_str_fi_intensitynew3_1_end<ppt){
            spe_ref_XCorr_yf_decoy_end[ncor_yf_decoy_end]<-0
            spe_yf_intensity_decoy_end[ncor_yf_decoy_end]<-0
          }else{
            cor_yf1_end<-vector()
            sum_str_fi_intensitynew1_end<-vector()
            for(npyf1_end in 1:nfimass){
              str_fi_intensitynew_end<-strsplit(str_fi_intensity_decoy_end[npyf1_end],";")
              str_fi_intensitynew1_end<-as.numeric(str_fi_intensitynew_end[[1]])
              nnozero_str_fi_intensitynew1_end<-sum(str_fi_intensitynew1_end!=0)
              if(nnozero_str_fi_intensitynew1_end<ppt){
                cor_yf1_end[npyf1_end]<-0
                sum_str_fi_intensitynew1_end[npyf1_end]<-0
              }else{
                cor_yf1_end[npyf1_end]<-sum(str_fi_intensitynew1_end*str_fi_intensitynew3_end)/sqrt(sum(str_fi_intensitynew1_end^2)*sum(str_fi_intensitynew3_end^2))
                sum_str_fi_intensitynew1_end[npyf1_end]<-log2(sum(str_fi_intensitynew1_end))
              }
            }
            if(all(cor_yf1_end==0)){
              spe_ref_XCorr_yf_decoy_end[ncor_yf_decoy_end]<-0
              spe_yf_intensity_decoy_end[ncor_yf_decoy_end]<-0
            }else{
              cor_yf_end_weight<-sum(cor_yf1_end)#*sum_str_fi_intensitynew1_end)/sum(sum_str_fi_intensitynew1_end)
              #cor_yf_noise_end<-str_fi_intensitynew3_end/str_fi_noisenew3_end
              #ncor_yf_noise_end<-sum(cor_yf_noise_end>0)
              #cor_yf_noise_index_end<-which(cor_yf_noise_end>0)
              #cor_yf_noise_sum_end<-sum(cor_yf_noise_end[cor_yf_noise_index_end])
              if(cor_yf_end_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yf_decoy_end[ncor_yf_decoy_end]<-cor_yf_end_weight
                nozero_str_fi_intensitynew3_end<-log2(sum(str_fi_intensitynew3_end))#*nnozero_str_fi_intensitynew3_end
                spe_yf_intensity_decoy_end[ncor_yf_decoy_end]<-nozero_str_fi_intensitynew3_end
              }else{
                spe_ref_XCorr_yf_decoy_end[ncor_yf_decoy_end]<-0
                spe_yf_intensity_decoy_end[ncor_yf_decoy_end]<-0
              }
              #*(1-exp(-(cor_yf_noise_sum_end-ncor_yf_noise_end)))
              
              
            }
          }
          ncor_yf_decoy_end<-ncor_yf_decoy_end+1
        }
        spe_ref_XCorr_yf_sum_decoy_end<-sum(spe_ref_XCorr_yf_decoy_end)
        spe_yf_intensity_sum_decoy_end<-sum(spe_yf_intensity_decoy_end)
        
        ####
        spe_yf_num_decoy_end<-sum(spe_ref_XCorr_yf_decoy_end!=0)
        
        if(max(spe_ref_XCorr_yf_decoy_end)>0){
          max_y_cor_index_end<-which.max(spe_ref_XCorr_yf_decoy_end)
          if(max_y_cor_index_end<=(nfi_yfmass_end-nfimass)/3){
            newmax_cor_by_decoy_end<-ylabel1[nby-nphk+(nfi_yfmass_end-nfimass)/3-max_y_cor_index_end+1]
          }
          else if(max_y_cor_index_end<=2*(nfi_yfmass_end-nfimass)/3 & max_y_cor_index_end>(nfi_yfmass_end-nfimass)/3){
            max_y_cor_index1_end<-max_y_cor_index_end-(nfi_yfmass_end-nfimass)/3
            newmax_cor_by_decoy_end<-ylabel2[nby-nphk+(nfi_yfmass_end-nfimass)/3-max_y_cor_index1_end+1]
          }
          else{
            max_y_cor_index1_end<-max_y_cor_index_end-2*(nfi_yfmass_end-nfimass)/3
            newmax_cor_by_decoy_end<-ylabel3[nby-nphk+(nfi_yfmass_end-nfimass)/3-max_y_cor_index1_end+1]
          }
        }else{
          newmax_cor_by_decoy_end<-"NA"
        }
        
        ####
        spe_ref_XCorr_sum_ph[phk]<-spe_ref_XCorr_bh_sum_end+spe_ref_XCorr_yf_sum_end
        spe_fi_num_ph[phk]<-spe_bh_num_end+spe_yf_num_end
        spe_fi_intensity_sum_ph[phk]<-spe_bh_intensity_sum_end+spe_yf_intensity_sum_end
        max_cor_by_ph[phk]<-paste(max_cor_by_ph_bh_end,newmax_cor_by_end,sep=" ")
        ####
        spe_ref_XCorr_sum_ph_decoy[phk]<-spe_ref_XCorr_bh_sum_decoy_end+spe_ref_XCorr_yf_sum_decoy_end
        spe_fi_num_ph_decoy[phk]<-spe_bh_num_decoy_end+spe_yf_num_decoy_end
        spe_fi_intensity_sum_ph_decoy[phk]<-spe_bh_intensity_sum_decoy_end+spe_yf_intensity_sum_decoy_end
        max_cor_by_ph_decoy[phk]<-paste(max_cor_by_ph_bh_decoy_end,newmax_cor_by_decoy_end,sep=" ")
      }
      else{
        pby_index1<-which(nphk==char_sty[[1]])
        phk_bmass_after11<-newbmass[nphk:(char_sty[[1]][pby_index1+1]-1)]
        phk_bmass_after12<-phk_bmass_after11-18.0152
        phk_bmass_after13<-phk_bmass_after11-17.0304
        phk_bmass_after1<-c(phk_bmass_after11,phk_bmass_after12,phk_bmass_after13)
        phk_bmass_ahead11<-newbmass[char_sty[[1]][pby_index1-1]:(nphk-1)]
        phk_bmass_ahead12<-phk_bmass_ahead11-18.0152
        phk_bmass_ahead13<-phk_bmass_ahead11-17.0304
        phk_bmass_ahead1<-c(phk_bmass_ahead11,phk_bmass_ahead12,phk_bmass_ahead13)
        phk_ymass_after11<-yyrev[(char_sty[[1]][pby_index1-1]+1):nphk]
        phk_ymass_after12<-phk_ymass_after11-18.0152
        phk_ymass_after13<-phk_ymass_after11-17.0304
        phk_ymass_after1<-c(phk_ymass_after11,phk_ymass_after12,phk_ymass_after13)
        phk_ymass_ahead11<-yyrev[(nphk+1):char_sty[[1]][pby_index1+1]]
        phk_ymass_ahead12<-phk_ymass_ahead11-18.0152
        phk_ymass_ahead13<-phk_ymass_ahead11-17.0304
        phk_ymass_ahead1<-c(phk_ymass_ahead11,phk_ymass_ahead12,phk_ymass_ahead13)
        fi_bfmass<-c(fimass,phk_bmass_after1)
        fi_bhmass<-c(fimass,phk_bmass_ahead1)
        fi_yfmass<-c(fimass,phk_ymass_after1)
        fi_yhmass<-c(fimass,phk_ymass_ahead1)
        nfi_bfmass<-length(fi_bfmass)
        nfi_bhmass<-length(fi_bhmass)
        nfi_yfmass<-length(fi_yfmass)
        nfi_yhmass<-length(fi_yhmass)
        bf_fimass_ph1_intensity<-vector()
        bh_fimass_ph1_intensity<-vector()
        yf_fimass_ph1_intensity<-vector()
        yh_fimass_ph1_intensity<-vector()
        bf_str_fi_intensity<-vector()
        bh_str_fi_intensity<-vector()
        yf_str_fi_intensity<-vector()
        yh_str_fi_intensity<-vector()
        bf_fimass_ph1_noise<-vector()
        bh_fimass_ph1_noise<-vector()
        yf_fimass_ph1_noise<-vector()
        yh_fimass_ph1_noise<-vector()
        bf_str_fi_noise<-vector()
        bh_str_fi_noise<-vector()
        yf_str_fi_noise<-vector()
        yh_str_fi_noise<-vector()
        ####
        phk_bmass_after_decoy11<-newbmass_decoy[nphk:(char_sty[[1]][pby_index1+1]-1)]
        phk_bmass_after_decoy12<-phk_bmass_after_decoy11-18.0152
        phk_bmass_after_decoy13<-phk_bmass_after_decoy11-17.0304
        phk_bmass_after_decoy1<-c(phk_bmass_after_decoy11,phk_bmass_after_decoy12,
                                  phk_bmass_after_decoy13)
        phk_bmass_ahead_decoy11<-newbmass_decoy[char_sty[[1]][pby_index1-1]:(nphk-1)]
        phk_bmass_ahead_decoy12<-phk_bmass_ahead_decoy11-18.0152
        phk_bmass_ahead_decoy13<-phk_bmass_ahead_decoy11-17.0304
        phk_bmass_ahead_decoy1<-c(phk_bmass_ahead_decoy11,phk_bmass_ahead_decoy12,
                                  phk_bmass_ahead_decoy13)
        phk_ymass_after_decoy11<-yyrev_decoy[(char_sty[[1]][pby_index1-1]+1):nphk]
        phk_ymass_after_decoy12<-phk_ymass_after_decoy11-18.0152
        phk_ymass_after_decoy13<-phk_ymass_after_decoy11-17.0304
        phk_ymass_after_decoy1<-c(phk_ymass_after_decoy11,phk_ymass_after_decoy12,
                                  phk_ymass_after_decoy13)
        phk_ymass_ahead_decoy11<-yyrev_decoy[(nphk+1):char_sty[[1]][pby_index1+1]]
        phk_ymass_ahead_decoy12<-phk_ymass_ahead_decoy11-18.0152
        phk_ymass_ahead_decoy13<-phk_ymass_ahead_decoy11-17.0304
        phk_ymass_ahead_decoy1<-c(phk_ymass_ahead_decoy11,phk_ymass_ahead_decoy12,
                                  phk_ymass_ahead_decoy13)
        fi_bfmass_decoy<-c(fimass_decoy,phk_bmass_after_decoy1)
        fi_bhmass_decoy<-c(fimass_decoy,phk_bmass_ahead_decoy1)
        fi_yfmass_decoy<-c(fimass_decoy,phk_ymass_after_decoy1)
        fi_yhmass_decoy<-c(fimass_decoy,phk_ymass_ahead_decoy1)
        nfi_bfmass_decoy<-length(fi_bfmass_decoy)
        nfi_bhmass_decoy<-length(fi_bhmass_decoy)
        nfi_yfmass_decoy<-length(fi_yfmass_decoy)
        nfi_yhmass_decoy<-length(fi_yhmass_decoy)
        bf_fimass_ph1_intensity_decoy<-vector()
        bh_fimass_ph1_intensity_decoy<-vector()
        yf_fimass_ph1_intensity_decoy<-vector()
        yh_fimass_ph1_intensity_decoy<-vector()
        bf_str_fi_intensity_decoy<-vector()
        bh_str_fi_intensity_decoy<-vector()
        yf_str_fi_intensity_decoy<-vector()
        yh_str_fi_intensity_decoy<-vector()
        bf_fimass_ph1_noise_decoy<-vector()
        bh_fimass_ph1_noise_decoy<-vector()
        yf_fimass_ph1_noise_decoy<-vector()
        yh_fimass_ph1_noise_decoy<-vector()
        bf_str_fi_noise_decoy<-vector()
        bh_str_fi_noise_decoy<-vector()
        yf_str_fi_noise_decoy<-vector()
        yh_str_fi_noise_decoy<-vector()
        ####
        
        for(bfsfi1 in 1:nfi_bfmass){
          for(nrt in 1:nnewrt){
            nrt1<-newrt[nrt]
            char_mz<-as.character(datapsm01$m.z[nrt1])
            str_mz<-strsplit(char_mz,";")
            newstr_mz<-as.numeric(str_mz[[1]])
            char_intensity<-as.character(datapsm01$Intensity[nrt1])
            str_intensity<-strsplit(char_intensity,";")
            newstr_intensity<-as.numeric(str_intensity[[1]])
            char_noise<-as.character(datapsm01$Noise[nrt1])
            str_noise<-strsplit(char_noise,";")
            newstr_noise<-as.numeric(str_noise[[1]])
            
            sub_fion1<-abs(newstr_mz-fi_bfmass[bfsfi1])
            if(min(sub_fion1)<=peakTolThreshold){
              sub_fion1_min<-which.min(sub_fion1)
              bf_fimass_ph1_intensity[nrt]<-newstr_intensity[sub_fion1_min]
              bf_fimass_ph1_noise[nrt]<-newstr_noise[sub_fion1_min]
            }else{
              bf_fimass_ph1_intensity[nrt]<-0
              bf_fimass_ph1_noise[nrt]<-0.0001
            }
            ####
            sub_fion1_decoy<-abs(newstr_mz-fi_bfmass_decoy[bfsfi1])
            if(min(sub_fion1_decoy)<=peakTolThreshold){
              sub_fion1_min_decoy<-which.min(sub_fion1_decoy)
              bf_fimass_ph1_intensity_decoy[nrt]<-newstr_intensity[sub_fion1_min_decoy]
              bf_fimass_ph1_noise_decoy[nrt]<-newstr_noise[sub_fion1_min_decoy]
            }else{
              bf_fimass_ph1_intensity_decoy[nrt]<-0
              bf_fimass_ph1_noise_decoy[nrt]<-0.0001
            }
          }
          newfimass_ph1_inten<-paste(bf_fimass_ph1_intensity,collapse=";")
          bf_str_fi_intensity[bfsfi1]<-newfimass_ph1_inten
          newfimass_ph1_noise<-paste(bf_fimass_ph1_noise,collapse=";")
          bf_str_fi_noise[bfsfi1]<-newfimass_ph1_noise
          ####
          newfimass_ph1_inten_decoy<-paste(bf_fimass_ph1_intensity_decoy,collapse=";")
          bf_str_fi_intensity_decoy[bfsfi1]<-newfimass_ph1_inten_decoy
          newfimass_ph1_noise_decoy<-paste(bf_fimass_ph1_noise_decoy,collapse=";")
          bf_str_fi_noise_decoy[bfsfi1]<-newfimass_ph1_noise_decoy
        }
        
        #####Target_b ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bf<-vector()
        spe_bf_intensity<-vector()
        ncor_bf<-1
        for(npbf2 in (nfimass+1):nfi_bfmass){
          bf_str_fi_intensitynew2<-strsplit(bf_str_fi_intensity[npbf2],";")
          bf_str_fi_intensitynew3<-as.numeric(bf_str_fi_intensitynew2[[1]])
          nnozero_bf_str_fi_intensitynew3_1<-sum(bf_str_fi_intensitynew3!=0)
          bf_str_fi_noisenew2<-strsplit(bf_str_fi_noise[npbf2],";")
          bf_str_fi_noisenew3<-as.numeric(bf_str_fi_noisenew2[[1]])
          if(nnozero_bf_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_bf[ncor_bf]<-0
            spe_bf_intensity[ncor_bf]<-0
          }else{
            cor_bf1<-vector()
            sum_bf_str_fi_intensitynew1<-vector()
            for(npbf1 in 1:nfimass){
              bf_str_fi_intensitynew<-strsplit(bf_str_fi_intensity[npbf1],";")
              bf_str_fi_intensitynew1<-as.numeric(bf_str_fi_intensitynew[[1]])
              nnozero_bf_str_fi_intensitynew1<-sum(bf_str_fi_intensitynew1!=0)
              if(nnozero_bf_str_fi_intensitynew1<ppt){
                cor_bf1[npbf1]<-0
                sum_bf_str_fi_intensitynew1[npbf1]<-0
              }else{
                cor_bf1[npbf1]<-sum(bf_str_fi_intensitynew1*bf_str_fi_intensitynew3)/sqrt(sum(bf_str_fi_intensitynew1^2)*sum(bf_str_fi_intensitynew3^2))
                sum_bf_str_fi_intensitynew1[npbf1]<-log2(sum(bf_str_fi_intensitynew1))
              }
            }
            if(all(cor_bf1==0)){
              spe_ref_XCorr_bf[ncor_bf]<-0
              spe_bf_intensity[ncor_bf]<-0
            }else{
              cor_bf_weight<-sum(cor_bf1)#*sum_bf_str_fi_intensitynew1)/sum(sum_bf_str_fi_intensitynew1)
              #cor_bf_noise<-bf_str_fi_intensitynew3/bf_str_fi_noisenew3
              #ncor_bf_noise<-sum(cor_bf_noise>0)
              #cor_bf_noise_index<-which(cor_bf_noise>0)
              #cor_bf_noise_sum<-sum(cor_bf_noise[cor_bf_noise_index])
              if(cor_bf_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bf[ncor_bf]<-cor_bf_weight
                nozero_bf_str_fi_intensitynew3<-log2(sum(bf_str_fi_intensitynew3))#*nnozero_bf_str_fi_intensitynew3
                spe_bf_intensity[ncor_bf]<-nozero_bf_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_bf[ncor_bf]<-0
                spe_bf_intensity[ncor_bf]<-0
              }
              #*(1-exp(-(cor_bf_noise_sum-ncor_bf_noise)))
              
              
            }
          }
          ncor_bf<-ncor_bf+1
        }
        spe_ref_XCorr_bf_sum<-sum(spe_ref_XCorr_bf)
        spe_bf_intensity_sum<-sum(spe_bf_intensity)
        
        ####
        spe_bf_num<-sum(spe_ref_XCorr_bf!=0)
        
        if(max(spe_ref_XCorr_bf)>0){
          max_bf_cor_index<-which.max(spe_ref_XCorr_bf)
          if(max_bf_cor_index<=(nfi_bfmass-nfimass)/3){
            max_cor_bf<-blabel1[nphk+max_bf_cor_index-1]
          }
          else if(max_bf_cor_index<=2*(nfi_bfmass-nfimass)/3 & max_bf_cor_index>(nfi_bfmass-nfimass)/3){
            max_bf_cor_index1<-max_bf_cor_index-(nfi_bfmass-nfimass)/3
            max_cor_bf<-blabel2[nphk+max_bf_cor_index1-1]
          }
          else{
            max_bf_cor_index2<-max_bf_cor_index-2*(nfi_bfmass-nfimass)/3
            max_cor_bf<-blabel3[nphk+max_bf_cor_index2-1]
          }
        }else{
          max_cor_bf<-"NA"
        }
        
        #####Decoy_b ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bf_decoy<-vector()
        spe_bf_intensity_decoy<-vector()
        ncor_bf_decoy<-1
        for(npbf2 in (nfimass+1):nfi_bfmass){
          bf_str_fi_intensitynew2<-strsplit(bf_str_fi_intensity_decoy[npbf2],";")
          bf_str_fi_intensitynew3<-as.numeric(bf_str_fi_intensitynew2[[1]])
          nnozero_bf_str_fi_intensitynew3_1<-sum(bf_str_fi_intensitynew3!=0)
          bf_str_fi_noisenew2<-strsplit(bf_str_fi_noise_decoy[npbf2],";")
          bf_str_fi_noisenew3<-as.numeric(bf_str_fi_noisenew2[[1]])
          if(nnozero_bf_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_bf_decoy[ncor_bf_decoy]<-0
            spe_bf_intensity_decoy[ncor_bf_decoy]<-0
          }else{
            cor_bf1<-vector()
            sum_bf_str_fi_intensitynew1<-vector()
            for(npbf1 in 1:nfimass){
              bf_str_fi_intensitynew<-strsplit(bf_str_fi_intensity_decoy[npbf1],";")
              bf_str_fi_intensitynew1<-as.numeric(bf_str_fi_intensitynew[[1]])
              nnozero_bf_str_fi_intensitynew1<-sum(bf_str_fi_intensitynew1!=0)
              if(nnozero_bf_str_fi_intensitynew1<ppt){
                cor_bf1[npbf1]<-0
                sum_bf_str_fi_intensitynew1[npbf1]<-0
              }else{
                cor_bf1[npbf1]<-sum(bf_str_fi_intensitynew1*bf_str_fi_intensitynew3)/sqrt(sum(bf_str_fi_intensitynew1^2)*sum(bf_str_fi_intensitynew3^2))
                sum_bf_str_fi_intensitynew1[npbf1]<-log2(sum(bf_str_fi_intensitynew1))
              }
            }
            if(all(cor_bf1==0)){
              spe_ref_XCorr_bf_decoy[ncor_bf_decoy]<-0
              spe_bf_intensity_decoy[ncor_bf_decoy]<-0
            }else{
              cor_bf_weight<-sum(cor_bf1)#*sum_bf_str_fi_intensitynew1)/sum(sum_bf_str_fi_intensitynew1)
              #cor_bf_noise<-bf_str_fi_intensitynew3/bf_str_fi_noisenew3
              #ncor_bf_noise<-sum(cor_bf_noise>0)
              #cor_bf_noise_index<-which(cor_bf_noise>0)
              #cor_bf_noise_sum<-sum(cor_bf_noise[cor_bf_noise_index])
              if( cor_bf_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bf_decoy[ncor_bf_decoy]<-cor_bf_weight
                nozero_bf_str_fi_intensitynew3<-log2(sum(bf_str_fi_intensitynew3))#*nnozero_bf_str_fi_intensitynew3
                spe_bf_intensity_decoy[ncor_bf_decoy]<-nozero_bf_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_bf_decoy[ncor_bf_decoy]<-0
                spe_bf_intensity_decoy[ncor_bf_decoy]<-0
              }
              #*(1-exp(-(cor_bf_noise_sum-ncor_bf_noise)))
              
            }
          }
          ncor_bf_decoy<-ncor_bf_decoy+1
        }
        spe_ref_XCorr_bf_sum_decoy<-sum(spe_ref_XCorr_bf_decoy)
        spe_bf_intensity_sum_decoy<-sum(spe_bf_intensity_decoy)
        
        ####
        spe_bf_num_decoy<-sum(spe_ref_XCorr_bf_decoy!=0)
        
        if(max(spe_ref_XCorr_bf_decoy)>0){
          max_bf_cor_index<-which.max(spe_ref_XCorr_bf_decoy)
          if(max_bf_cor_index<=(nfi_bfmass-nfimass)/3){
            max_cor_bf_decoy<-blabel1[nphk+max_bf_cor_index-1]
          }
          else if(max_bf_cor_index<=2*(nfi_bfmass-nfimass)/3 & max_bf_cor_index>(nfi_bfmass-nfimass)/3){
            max_bf_cor_index1<-max_bf_cor_index-(nfi_bfmass-nfimass)/3
            max_cor_bf_decoy<-blabel2[nphk+max_bf_cor_index1-1]
          }
          else{
            max_bf_cor_index2<-max_bf_cor_index-2*(nfi_bfmass-nfimass)/3
            max_cor_bf_decoy<-blabel3[nphk+max_bf_cor_index2-1]
          }
        }else{
          max_cor_bf_decoy<-"NA"
        }
        
        ####
        #
        for(bhsfi1 in 1:nfi_bhmass){
          for(nrt in 1:nnewrt){
            nrt1<-newrt[nrt]
            char_mz<-as.character(datapsm01$m.z[nrt1])
            str_mz<-strsplit(char_mz,";")
            newstr_mz<-as.numeric(str_mz[[1]])
            char_intensity<-as.character(datapsm01$Intensity[nrt1])
            str_intensity<-strsplit(char_intensity,";")
            newstr_intensity<-as.numeric(str_intensity[[1]])
            char_noise<-as.character(datapsm01$Noise[nrt1])
            str_noise<-strsplit(char_noise,";")
            newstr_noise<-as.numeric(str_noise[[1]])
            
            sub_fion1<-abs(newstr_mz-fi_bhmass[bhsfi1])
            if(min(sub_fion1)<=peakTolThreshold){
              sub_fion1_min<-which.min(sub_fion1)
              bh_fimass_ph1_intensity[nrt]<-newstr_intensity[sub_fion1_min]
              bh_fimass_ph1_noise[nrt]<-newstr_noise[sub_fion1_min]
            }else{
              bh_fimass_ph1_intensity[nrt]<-0
              bh_fimass_ph1_noise[nrt]<-0.0001
            }
            ####
            sub_fion1_decoy<-abs(newstr_mz-fi_bhmass_decoy[bhsfi1])
            if(min(sub_fion1_decoy)<=peakTolThreshold){
              sub_fion1_min_decoy<-which.min(sub_fion1_decoy)
              bh_fimass_ph1_intensity_decoy[nrt]<-newstr_intensity[sub_fion1_min_decoy]
              bh_fimass_ph1_noise_decoy[nrt]<-newstr_noise[sub_fion1_min_decoy]
            }else{
              bh_fimass_ph1_intensity_decoy[nrt]<-0
              bh_fimass_ph1_noise_decoy[nrt]<-0.0001
            }
          }
          newfimass_ph1_inten<-paste(bh_fimass_ph1_intensity,collapse=";")
          bh_str_fi_intensity[bhsfi1]<-newfimass_ph1_inten
          newfimass_ph1_noise<-paste(bh_fimass_ph1_noise,collapse=";")
          bh_str_fi_noise[bhsfi1]<-newfimass_ph1_noise
          ####
          newfimass_ph1_inten_decoy<-paste(bh_fimass_ph1_intensity_decoy,collapse=";")
          bh_str_fi_intensity_decoy[bhsfi1]<-newfimass_ph1_inten_decoy
          newfimass_ph1_noise_decoy<-paste(bh_fimass_ph1_noise_decoy,collapse=";")
          bh_str_fi_noise_decoy[bhsfi1]<-newfimass_ph1_noise_decoy
        }
        
        #####Target_b ion ahead
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bh<-vector()
        spe_bh_intensity<-vector()
        ncor_bh<-1
        for(npbh2 in (nfimass+1):nfi_bhmass){
          bh_str_fi_intensitynew2<-strsplit(bh_str_fi_intensity[npbh2],";")
          bh_str_fi_intensitynew3<-as.numeric(bh_str_fi_intensitynew2[[1]])
          nnozero_bh_str_fi_intensitynew3_1<-sum(bh_str_fi_intensitynew3!=0)
          bh_str_fi_noisenew2<-strsplit(bh_str_fi_noise[npbh2],";")
          bh_str_fi_noisenew3<-as.numeric(bh_str_fi_noisenew2[[1]])
          if(nnozero_bh_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_bh[ncor_bh]<-0
            spe_bh_intensity[ncor_bh]<-0
          }else{
            cor_bh1<-vector()
            sum_bh_str_fi_intensitynew1<-vector()
            for(npbh1 in 1:nfimass){
              bh_str_fi_intensitynew<-strsplit(bh_str_fi_intensity[npbh1],";")
              bh_str_fi_intensitynew1<-as.numeric(bh_str_fi_intensitynew[[1]])
              nnozero_bh_str_fi_intensitynew1<-sum(bh_str_fi_intensitynew1!=0)
              if(nnozero_bh_str_fi_intensitynew1<ppt){
                cor_bh1[npbh1]<-0
                sum_bh_str_fi_intensitynew1[npbh1]<-0
              }else{
                cor_bh1[npbh1]<-sum(bh_str_fi_intensitynew1*bh_str_fi_intensitynew3)/sqrt(sum(bh_str_fi_intensitynew1^2)*sum(bh_str_fi_intensitynew3^2))
                sum_bh_str_fi_intensitynew1[npbh1]<-log2(sum(bh_str_fi_intensitynew1))
              }
            }
            if(all(cor_bh1==0)){
              spe_ref_XCorr_bh[ncor_bh]<-0
              spe_bh_intensity[ncor_bh]<-0
            }else{
              cor_bh_weight<-sum(cor_bh1)#*sum_bh_str_fi_intensitynew1)/sum(sum_bh_str_fi_intensitynew1)
              #cor_bh_noise<-bh_str_fi_intensitynew3/bh_str_fi_noisenew3
              #ncor_bh_noise<-sum(cor_bh_noise>0)
              #cor_bh_noise_index<-which(cor_bh_noise>0)
              #cor_bh_noise_sum<-sum(cor_bh_noise[cor_bh_noise_index])
              if(cor_bh_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bh[ncor_bh]<-cor_bh_weight
                nozero_bh_str_fi_intensitynew3<-log2(sum(bh_str_fi_intensitynew3))#*nnozero_bh_str_fi_intensitynew3
                spe_bh_intensity[ncor_bh]<-nozero_bh_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_bh[ncor_bh]<-0
                spe_bh_intensity[ncor_bh]<-0
              }
              #*(1-exp(-(cor_bh_noise_sum-ncor_bh_noise)))
              
            }
          }
          ncor_bh<-ncor_bh+1
        }
        spe_ref_XCorr_bh_sum<-sum(spe_ref_XCorr_bh)
        spe_bh_intensity_sum<-sum(spe_bh_intensity)
        
        ####
        spe_bh_num<-sum(spe_ref_XCorr_bh!=0)
        
        if(max(spe_ref_XCorr_bh)>0){
          max_bh_cor_index<-which.max(spe_ref_XCorr_bh)
          if(max_bh_cor_index<=(nfi_bhmass-nfimass)/3){
            max_cor_bh<-blabel1[nphk-(nfi_bhmass-nfimass)/3+max_bh_cor_index-1]
          }
          else if(max_bh_cor_index<=2*(nfi_bhmass-nfimass)/3 & max_bh_cor_index>(nfi_bhmass-nfimass)/3){
            max_bh_cor_index1<-max_bh_cor_index-(nfi_bhmass-nfimass)/3
            max_cor_bh<-blabel2[nphk-(nfi_bhmass-nfimass)/3+max_bh_cor_index1-1]
          }
          else{
            max_bh_cor_index2<-max_bh_cor_index-2*(nfi_bhmass-nfimass)/3
            max_cor_bh<-blabel3[nphk-(nfi_bhmass-nfimass)/3+max_bh_cor_index2-1]
          }
        }else{
          max_cor_bh<-"NA"
        }
        
        #####Decoy_b ion ahead
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_bh_decoy<-vector()
        spe_bh_intensity_decoy<-vector()
        ncor_bh_decoy<-1
        for(npbh2 in (nfimass+1):nfi_bhmass){
          bh_str_fi_intensitynew2<-strsplit(bh_str_fi_intensity_decoy[npbh2],";")
          bh_str_fi_intensitynew3<-as.numeric(bh_str_fi_intensitynew2[[1]])
          nnozero_bh_str_fi_intensitynew3_1<-sum(bh_str_fi_intensitynew3!=0)
          bh_str_fi_noisenew2<-strsplit(bh_str_fi_noise_decoy[npbh2],";")
          bh_str_fi_noisenew3<-as.numeric(bh_str_fi_noisenew2[[1]])
          if(nnozero_bh_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_bh_decoy[ncor_bh_decoy]<-0
            spe_bh_intensity_decoy[ncor_bh_decoy]<-0
          }else{
            cor_bh1<-vector()
            sum_bh_str_fi_intensitynew1<-vector()
            for(npbh1 in 1:nfimass){
              bh_str_fi_intensitynew<-strsplit(bh_str_fi_intensity_decoy[npbh1],";")
              bh_str_fi_intensitynew1<-as.numeric(bh_str_fi_intensitynew[[1]])
              nnozero_bh_str_fi_intensitynew1<-sum(bh_str_fi_intensitynew1!=0)
              if(nnozero_bh_str_fi_intensitynew1<ppt){
                cor_bh1[npbh1]<-0
                sum_bh_str_fi_intensitynew1[npbh1]<-0
              }else{
                cor_bh1[npbh1]<-sum(bh_str_fi_intensitynew1*bh_str_fi_intensitynew3)/sqrt(sum(bh_str_fi_intensitynew1^2)*sum(bh_str_fi_intensitynew3^2))
                sum_bh_str_fi_intensitynew1[npbh1]<-log2(sum(bh_str_fi_intensitynew1))
              }
            }
            if(all(cor_bh1==0)){
              spe_ref_XCorr_bh_decoy[ncor_bh_decoy]<-0
              spe_bh_intensity_decoy[ncor_bh_decoy]<-0
            }else{
              cor_bh_weight<-sum(cor_bh1)#*sum_bh_str_fi_intensitynew1)/sum(sum_bh_str_fi_intensitynew1)
              #cor_bh_noise<-bh_str_fi_intensitynew3/bh_str_fi_noisenew3
              #ncor_bh_noise<-sum(cor_bh_noise>0)
              #cor_bh_noise_index<-which(cor_bh_noise>0)
              #cor_bh_noise_sum<-sum(cor_bh_noise[cor_bh_noise_index])
              if(cor_bh_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_bh_decoy[ncor_bh_decoy]<-cor_bh_weight
                nozero_bh_str_fi_intensitynew3<-log2(sum(bh_str_fi_intensitynew3))#*nnozero_bh_str_fi_intensitynew3
                spe_bh_intensity_decoy[ncor_bh_decoy]<-nozero_bh_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_bh_decoy[ncor_bh_decoy]<-0
                spe_bh_intensity_decoy[ncor_bh_decoy]<-0
              }
              #*(1-exp(-(cor_bh_noise_sum-ncor_bh_noise)))nnozero_bh_str_fi_intensitynew3<-sum(bh_str_fi_intensitynew3!=0)
              
            }
          }
          ncor_bh_decoy<-ncor_bh_decoy+1
        }
        spe_ref_XCorr_bh_sum_decoy<-sum(spe_ref_XCorr_bh_decoy)
        spe_bh_intensity_sum_decoy<-sum(spe_bh_intensity_decoy)
        
        ####
        spe_bh_num_decoy<-sum(spe_ref_XCorr_bh_decoy!=0)
        
        if(max(spe_ref_XCorr_bh_decoy)>0){
          max_bh_cor_index<-which.max(spe_ref_XCorr_bh_decoy)
          if(max_bh_cor_index<=(nfi_bhmass-nfimass)/3){
            max_cor_bh_decoy<-blabel1[nphk-(nfi_bhmass-nfimass)/3+max_bh_cor_index-1]
          }
          else if(max_bh_cor_index<=2*(nfi_bhmass-nfimass)/3 & max_bh_cor_index>(nfi_bhmass-nfimass)/3){
            max_bh_cor_index1<-max_bh_cor_index-(nfi_bhmass-nfimass)/3
            max_cor_bh_decoy<-blabel2[nphk-(nfi_bhmass-nfimass)/3+max_bh_cor_index1-1]
          }
          else{
            max_bh_cor_index2<-max_bh_cor_index-2*(nfi_bhmass-nfimass)/3
            max_cor_bh_decoy<-blabel3[nphk-(nfi_bhmass-nfimass)/3+max_bh_cor_index2-1]
          }
        }else{
          max_cor_bh_decoy<-"NA"
        }
        
        ####
        #
        for(yfsfi1 in 1:nfi_yfmass){
          for(nrt in 1:nnewrt){
            nrt1<-newrt[nrt]
            char_mz<-as.character(datapsm01$m.z[nrt1])
            str_mz<-strsplit(char_mz,";")
            newstr_mz<-as.numeric(str_mz[[1]])
            char_intensity<-as.character(datapsm01$Intensity[nrt1])
            str_intensity<-strsplit(char_intensity,";")
            newstr_intensity<-as.numeric(str_intensity[[1]])
            char_noise<-as.character(datapsm01$Noise[nrt1])
            str_noise<-strsplit(char_noise,";")
            newstr_noise<-as.numeric(str_noise[[1]])
            
            sub_fion1<-abs(newstr_mz-fi_yfmass[yfsfi1])
            if(min(sub_fion1)<=peakTolThreshold){
              sub_fion1_min<-which.min(sub_fion1)
              yf_fimass_ph1_intensity[nrt]<-newstr_intensity[sub_fion1_min]
              yf_fimass_ph1_noise[nrt]<-newstr_noise[sub_fion1_min]
            }
            else{
              yf_fimass_ph1_intensity[nrt]<-0
              yf_fimass_ph1_noise[nrt]<-0.0001
            }
            ####
            sub_fion1_decoy<-abs(newstr_mz-fi_yfmass_decoy[yfsfi1])
            if(min(sub_fion1_decoy)<=peakTolThreshold){
              sub_fion1_min_decoy<-which.min(sub_fion1_decoy)
              yf_fimass_ph1_intensity_decoy[nrt]<-newstr_intensity[sub_fion1_min_decoy]
              yf_fimass_ph1_noise_decoy[nrt]<-newstr_noise[sub_fion1_min_decoy]
            }else{
              yf_fimass_ph1_intensity_decoy[nrt]<-0
              yf_fimass_ph1_noise_decoy[nrt]<-0.0001
            }
          }
          newfimass_ph1_inten<-paste(yf_fimass_ph1_intensity,collapse=";")
          yf_str_fi_intensity[yfsfi1]<-newfimass_ph1_inten
          newfimass_ph1_noise<-paste(yf_fimass_ph1_noise,collapse=";")
          yf_str_fi_noise[yfsfi1]<-newfimass_ph1_noise
          ####
          newfimass_ph1_inten_decoy<-paste(yf_fimass_ph1_intensity_decoy,collapse=";")
          yf_str_fi_intensity_decoy[yfsfi1]<-newfimass_ph1_inten_decoy
          newfimass_ph1_noise_decoy<-paste(yf_fimass_ph1_noise_decoy,collapse=";")
          yf_str_fi_noise_decoy[yfsfi1]<-newfimass_ph1_noise_decoy
        }
        
        #####Target_y ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yf<-vector()
        spe_yf_intensity<-vector()
        ncor_yf<-1
        for(npyf2 in (nfimass+1):nfi_yfmass){
          yf_str_fi_intensitynew2<-strsplit(yf_str_fi_intensity[npyf2],";")
          yf_str_fi_intensitynew3<-as.numeric(yf_str_fi_intensitynew2[[1]])
          nnozero_yf_str_fi_intensitynew3_1<-sum(yf_str_fi_intensitynew3!=0)
          yf_str_fi_noisenew2<-strsplit(yf_str_fi_noise[npyf2],";")
          yf_str_fi_noisenew3<-as.numeric(yf_str_fi_noisenew2[[1]])
          if(nnozero_yf_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_yf[ncor_yf]<-0
            spe_yf_intensity[ncor_yf]<-0
          }else{
            cor_yf1<-vector()
            sum_yf_str_fi_intensitynew1<-vector()
            for(npyf1 in 1:nfimass){
              yf_str_fi_intensitynew<-strsplit(yf_str_fi_intensity[npyf1],";")
              yf_str_fi_intensitynew1<-as.numeric(yf_str_fi_intensitynew[[1]])
              nnozero_yf_str_fi_intensitynew1<-sum(yf_str_fi_intensitynew1!=0)
              if(nnozero_yf_str_fi_intensitynew1<ppt){
                cor_yf1[npyf1]<-0
                sum_yf_str_fi_intensitynew1[npyf1]<-0
              }else{
                cor_yf1[npyf1]<-sum(yf_str_fi_intensitynew1*yf_str_fi_intensitynew3)/sqrt(sum(yf_str_fi_intensitynew1^2)*sum(yf_str_fi_intensitynew3^2))
                sum_yf_str_fi_intensitynew1[npyf1]<-log2(sum(yf_str_fi_intensitynew1))
              }
            }
            if(all(cor_yf1==0)){
              spe_ref_XCorr_yf[ncor_yf]<-0
              spe_yf_intensity[ncor_yf]<-0
            }else{
              cor_yf_weight<-sum(cor_yf1)#*sum_yf_str_fi_intensitynew1)/sum(sum_yf_str_fi_intensitynew1)
              #cor_yf_noise<-yf_str_fi_intensitynew3/yf_str_fi_noisenew3
              #ncor_yf_noise<-sum(cor_yf_noise>0)
              #cor_yf_noise_index<-which(cor_yf_noise>0)
              #cor_yf_noise_sum<-sum(cor_yf_noise[cor_yf_noise_index])
              if(cor_yf_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yf[ncor_yf]<-cor_yf_weight
                nozero_yf_str_fi_intensitynew3<-log2(sum(yf_str_fi_intensitynew3))#*nnozero_yf_str_fi_intensitynew3
                spe_yf_intensity[ncor_yf]<-nozero_yf_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_yf[ncor_yf]<-0
                spe_yf_intensity[ncor_yf]<-0
              }
              #*(1-exp(-(cor_yf_noise_sum-ncor_yf_noise)))nnozero_yf_str_fi_intensitynew3<-sum(yf_str_fi_intensitynew3!=0)
              
            }
          }
          ncor_yf<-ncor_yf+1
        }
        spe_ref_XCorr_yf_sum<-sum(spe_ref_XCorr_yf)
        spe_yf_intensity_sum<-sum(spe_yf_intensity)
        
        ####
        spe_yf_num<-sum(spe_ref_XCorr_yf!=0)
        
        if(max(spe_ref_XCorr_yf)>0){
          max_yf_cor_index<-which.max(spe_ref_XCorr_yf)
          if(max_yf_cor_index<=(nfi_yfmass-nfimass)/3){
            max_cor_yf<-ylabel1[nby-nphk+(nfi_yfmass-nfimass)/3-max_yf_cor_index+1]
          }
          else if(max_yf_cor_index<=2*(nfi_yfmass-nfimass)/3 & max_yf_cor_index>(nfi_yfmass-nfimass)/3){
            max_yf_cor_index1<-max_yf_cor_index-(nfi_yfmass-nfimass)/3
            max_cor_yf<-ylabel2[nby-nphk+(nfi_yfmass-nfimass)/3-max_yf_cor_index1+1]
          }
          else{
            max_yf_cor_index2<-max_yf_cor_index-2*(nfi_yfmass-nfimass)/3
            max_cor_yf<-ylabel3[nby-nphk+(nfi_yfmass-nfimass)/3-max_yf_cor_index2+1]
          }
        }else{
          max_cor_yf<-"NA"
        }
        
        #####Decoy_y ion after
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yf_decoy<-vector()
        spe_yf_intensity_decoy<-vector()
        ncor_yf_decoy<-1
        for(npyf2 in (nfimass+1):nfi_yfmass){
          yf_str_fi_intensitynew2<-strsplit(yf_str_fi_intensity_decoy[npyf2],";")
          yf_str_fi_intensitynew3<-as.numeric(yf_str_fi_intensitynew2[[1]])
          nnozero_yf_str_fi_intensitynew3_1<-sum(yf_str_fi_intensitynew3!=0)
          yf_str_fi_noisenew2<-strsplit(yf_str_fi_noise_decoy[npyf2],";")
          yf_str_fi_noisenew3<-as.numeric(yf_str_fi_noisenew2[[1]])
          if(nnozero_yf_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_yf_decoy[ncor_yf_decoy]<-0
            spe_yf_intensity_decoy[ncor_yf_decoy]<-0
          }else{
            cor_yf1<-vector()
            sum_yf_str_fi_intensitynew1<-vector()
            for(npyf1 in 1:nfimass){
              yf_str_fi_intensitynew<-strsplit(yf_str_fi_intensity_decoy[npyf1],";")
              yf_str_fi_intensitynew1<-as.numeric(yf_str_fi_intensitynew[[1]])
              nnozero_yf_str_fi_intensitynew1<-sum(yf_str_fi_intensitynew1!=0)
              if(nnozero_yf_str_fi_intensitynew1<ppt){
                cor_yf1[npyf1]<-0
                sum_yf_str_fi_intensitynew1[npyf1]<-0
              }else{
                cor_yf1[npyf1]<-sum(yf_str_fi_intensitynew1*yf_str_fi_intensitynew3)/sqrt(sum(yf_str_fi_intensitynew1^2)*sum(yf_str_fi_intensitynew3^2))
                sum_yf_str_fi_intensitynew1[npyf1]<-log2(sum(yf_str_fi_intensitynew1))
              }
            }
            if(all(cor_yf1==0)){
              spe_ref_XCorr_yf_decoy[ncor_yf_decoy]<-0
              spe_yf_intensity_decoy[ncor_yf_decoy]<-0
            }else{
              cor_yf_weight<-sum(cor_yf1)#*sum_yf_str_fi_intensitynew1)/sum(sum_yf_str_fi_intensitynew1)
              #cor_yf_noise<-yf_str_fi_intensitynew3/yf_str_fi_noisenew3
              #ncor_yf_noise<-sum(cor_yf_noise>0)
              #cor_yf_noise_index<-which(cor_yf_noise>0)
              #cor_yf_noise_sum<-sum(cor_yf_noise[cor_yf_noise_index])
              if(cor_yf_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yf_decoy[ncor_yf_decoy]<-cor_yf_weight
                nozero_yf_str_fi_intensitynew3<-log2(sum(yf_str_fi_intensitynew3))#*nnozero_yf_str_fi_intensitynew3
                spe_yf_intensity_decoy[ncor_yf_decoy]<-nozero_yf_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_yf_decoy[ncor_yf_decoy]<-0
                spe_yf_intensity_decoy[ncor_yf_decoy]<-0
              }
              #*(1-exp(-(cor_yf_noise_sum-ncor_yf_noise)))
              
            }
          }
          ncor_yf_decoy<-ncor_yf_decoy+1
        }
        spe_ref_XCorr_yf_sum_decoy<-sum(spe_ref_XCorr_yf_decoy)
        spe_yf_intensity_sum_decoy<-sum(spe_yf_intensity_decoy)
        
        ####
        spe_yf_num_decoy<-sum(spe_ref_XCorr_yf_decoy!=0)
        
        if(max(spe_ref_XCorr_yf_decoy)>0){
          max_yf_cor_index<-which.max(spe_ref_XCorr_yf_decoy)
          if(max_yf_cor_index<=(nfi_yfmass-nfimass)/3){
            max_cor_yf_decoy<-ylabel1[nby-nphk+(nfi_yfmass-nfimass)/3-max_yf_cor_index+1]
          }
          else if(max_yf_cor_index<=2*(nfi_yfmass-nfimass)/3 & max_yf_cor_index>(nfi_yfmass-nfimass)/3){
            max_yf_cor_index1<-max_yf_cor_index-(nfi_yfmass-nfimass)/3
            max_cor_yf_decoy<-ylabel2[nby-nphk+(nfi_yfmass-nfimass)/3-max_yf_cor_index1+1]
          }
          else{
            max_yf_cor_index2<-max_yf_cor_index-2*(nfi_yfmass-nfimass)/3
            max_cor_yf_decoy<-ylabel3[nby-nphk+(nfi_yfmass-nfimass)/3-max_yf_cor_index2+1]
          }
        }else{
          max_cor_yf_decoy<-"NA"
        }
        
        ####
        #
        for(yhsfi1 in 1:nfi_yhmass){
          for(nrt in 1:nnewrt){
            nrt1<-newrt[nrt]
            char_mz<-as.character(datapsm01$m.z[nrt1])
            str_mz<-strsplit(char_mz,";")
            newstr_mz<-as.numeric(str_mz[[1]])
            char_intensity<-as.character(datapsm01$Intensity[nrt1])
            str_intensity<-strsplit(char_intensity,";")
            newstr_intensity<-as.numeric(str_intensity[[1]])
            char_noise<-as.character(datapsm01$Noise[nrt1])
            str_noise<-strsplit(char_noise,";")
            newstr_noise<-as.numeric(str_noise[[1]])
            
            sub_fion1<-abs(newstr_mz-fi_yhmass[yhsfi1])
            if(min(sub_fion1)<=peakTolThreshold){
              sub_fion1_min<-which.min(sub_fion1)
              yh_fimass_ph1_intensity[nrt]<-newstr_intensity[sub_fion1_min]
              yh_fimass_ph1_noise[nrt]<-newstr_noise[sub_fion1_min]
            }
            else{
              yh_fimass_ph1_intensity[nrt]<-0
              yh_fimass_ph1_noise[nrt]<-0.0001
            }
            ####
            sub_fion1_decoy<-abs(newstr_mz-fi_yhmass_decoy[yhsfi1])
            if(min(sub_fion1_decoy)<=peakTolThreshold){
              sub_fion1_min_decoy<-which.min(sub_fion1_decoy)
              yh_fimass_ph1_intensity_decoy[nrt]<-newstr_intensity[sub_fion1_min_decoy]
              yh_fimass_ph1_noise_decoy[nrt]<-newstr_noise[sub_fion1_min_decoy]
            }else{
              yh_fimass_ph1_intensity_decoy[nrt]<-0
              yh_fimass_ph1_noise_decoy[nrt]<-0.0001
            }
          }
          newfimass_ph1_inten<-paste(yh_fimass_ph1_intensity,collapse=";")
          yh_str_fi_intensity[yhsfi1]<-newfimass_ph1_inten
          newfimass_ph1_noise<-paste(yh_fimass_ph1_noise,collapse=";")
          yh_str_fi_noise[yhsfi1]<-newfimass_ph1_noise
          ####
          newfimass_ph1_inten_decoy<-paste(yh_fimass_ph1_intensity_decoy,collapse=";")
          yh_str_fi_intensity_decoy[yhsfi1]<-newfimass_ph1_inten_decoy
          newfimass_ph1_noise_decoy<-paste(yh_fimass_ph1_noise_decoy,collapse=";")
          yh_str_fi_noise_decoy[yhsfi1]<-newfimass_ph1_noise_decoy
        }
        
        #####Target_y ion ahead
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yh<-vector()
        spe_yh_intensity<-vector()
        ncor_yh<-1
        for(npyh2 in (nfimass+1):nfi_yhmass){
          yh_str_fi_intensitynew2<-strsplit(yh_str_fi_intensity[npyh2],";")
          yh_str_fi_intensitynew3<-as.numeric(yh_str_fi_intensitynew2[[1]])
          nnozero_yh_str_fi_intensitynew3_1<-sum(yh_str_fi_intensitynew3!=0)
          yh_str_fi_noisenew2<-strsplit(yh_str_fi_noise[npyh2],";")
          yh_str_fi_noisenew3<-as.numeric(yh_str_fi_noisenew2[[1]])
          if(nnozero_yh_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_yh[ncor_yh]<-0
            spe_yh_intensity[ncor_yh]<-0
          }else{
            cor_yh1<-vector()
            sum_yh_str_fi_intensitynew1<-vector()
            for(npyh1 in 1:nfimass){
              yh_str_fi_intensitynew<-strsplit(yh_str_fi_intensity[npyh1],";")
              yh_str_fi_intensitynew1<-as.numeric(yh_str_fi_intensitynew[[1]])
              nnozero_yh_str_fi_intensitynew1<-sum(yh_str_fi_intensitynew1!=0)
              if(nnozero_yh_str_fi_intensitynew1<ppt){
                cor_yh1[npyh1]<-0
                sum_yh_str_fi_intensitynew1[npyh1]<-0
              }else{
                cor_yh1[npyh1]<-sum(yh_str_fi_intensitynew1*yh_str_fi_intensitynew3)/sqrt(sum(yh_str_fi_intensitynew1^2)*sum(yh_str_fi_intensitynew3^2))
                sum_yh_str_fi_intensitynew1[npyh1]<-log2(sum(yh_str_fi_intensitynew1))
              }
            }
            if(all(cor_yh1==0)){
              spe_ref_XCorr_yh[ncor_yh]<-0
              spe_yh_intensity[ncor_yh]<-0
            }else{
              cor_yh_weight<-sum(cor_yh1)#*sum_yh_str_fi_intensitynew1)/sum(sum_yh_str_fi_intensitynew1)
              #cor_yh_noise<-yh_str_fi_intensitynew3/yh_str_fi_noisenew3
              #ncor_yh_noise<-sum(cor_yh_noise>0)
              #cor_yh_noise_index<-which(cor_yh_noise>0)
              #cor_yh_noise_sum<-sum(cor_yh_noise[cor_yh_noise_index])
              if(cor_yh_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yh[ncor_yh]<-cor_yh_weight
                nozero_yh_str_fi_intensitynew3<-log2(sum(yh_str_fi_intensitynew3))#*nnozero_yh_str_fi_intensitynew3
                spe_yh_intensity[ncor_yh]<-nozero_yh_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_yh[ncor_yh]<-0
                spe_yh_intensity[ncor_yh]<-0
              }
              #*(1-exp(-(cor_yh_noise_sum-ncor_yh_noise)))
              
            }
          }
          ncor_yh<-ncor_yh+1
        }
        spe_ref_XCorr_yh_sum<-sum(spe_ref_XCorr_yh)
        spe_yh_intensity_sum<-sum(spe_yh_intensity)
        
        ####
        spe_yh_num<-sum(spe_ref_XCorr_yh!=0)
        
        if(max(spe_ref_XCorr_yh)>0){
          max_yh_cor_index<-which.max(spe_ref_XCorr_yh)
          if(max_yh_cor_index<=(nfi_yhmass-nfimass)/3){
            max_cor_yh<-ylabel1[nby-nphk-max_yh_cor_index+1]
          }
          else if(max_yh_cor_index<=2*(nfi_yhmass-nfimass)/3 & max_yh_cor_index>(nfi_yhmass-nfimass)/3){
            max_yh_cor_index1<-max_yh_cor_index-(nfi_yhmass-nfimass)/3
            max_cor_yh<-ylabel2[nby-nphk-max_yh_cor_index1+1]
          }
          else{
            max_yh_cor_index2<-max_yh_cor_index-2*(nfi_yhmass-nfimass)/3
            max_cor_yh<-ylabel3[nby-nphk-max_yh_cor_index2+1]
          }
        }else{
          max_cor_yh<-"NA"
        }
        
        #####Decoy_y ion ahead
        ####Calculate the referencial ions autocorrelation coefficient
        
        
        ####Calculate the specific ions autocorrelation coefficient
        
        
        ####Calculate the cross-correlation of specific and referencial ions
        spe_ref_XCorr_yh_decoy<-vector()
        spe_yh_intensity_decoy<-vector()
        ncor_yh_decoy<-1
        for(npyh2 in (nfimass+1):nfi_yhmass){
          yh_str_fi_intensitynew2<-strsplit(yh_str_fi_intensity_decoy[npyh2],";")
          yh_str_fi_intensitynew3<-as.numeric(yh_str_fi_intensitynew2[[1]])
          nnozero_yh_str_fi_intensitynew3_1<-sum(yh_str_fi_intensitynew3!=0)
          yh_str_fi_noisenew2<-strsplit(yh_str_fi_noise_decoy[npyh2],";")
          yh_str_fi_noisenew3<-as.numeric(yh_str_fi_noisenew2[[1]])
          if(nnozero_yh_str_fi_intensitynew3_1<ppt){
            spe_ref_XCorr_yh_decoy[ncor_yh_decoy]<-0
            spe_yh_intensity_decoy[ncor_yh_decoy]<-0
          }else{
            cor_yh1<-vector()
            sum_yh_str_fi_intensitynew1<-vector()
            for(npyh1 in 1:nfimass){
              yh_str_fi_intensitynew<-strsplit(yh_str_fi_intensity_decoy[npyh1],";")
              yh_str_fi_intensitynew1<-as.numeric(yh_str_fi_intensitynew[[1]])
              nnozero_yh_str_fi_intensitynew1<-sum(yh_str_fi_intensitynew1!=0)
              if(nnozero_yh_str_fi_intensitynew1<ppt){
                cor_yh1[npyh1]<-0
                sum_yh_str_fi_intensitynew1[npyh1]<-0
              }else{
                cor_yh1[npyh1]<-sum(yh_str_fi_intensitynew1*yh_str_fi_intensitynew3)/sqrt(sum(yh_str_fi_intensitynew1^2)*sum(yh_str_fi_intensitynew3^2))
                sum_yh_str_fi_intensitynew1[npyh1]<-log2(sum(yh_str_fi_intensitynew1))
              }
            }
            if(all(cor_yh1==0)){
              spe_ref_XCorr_yh_decoy[ncor_yh_decoy]<-0
              spe_yh_intensity_decoy[ncor_yh_decoy]<-0
            }else{
              cor_yh_weight<-sum(cor_yh1)#*sum_yh_str_fi_intensitynew1)/sum(sum_yh_str_fi_intensitynew1)
              #cor_yh_noise<-yh_str_fi_intensitynew3/yh_str_fi_noisenew3
              #ncor_yh_noise<-sum(cor_yh_noise>0)
              #cor_yh_noise_index<-which(cor_yh_noise>0)
              #cor_yh_noise_sum<-sum(cor_yh_noise[cor_yh_noise_index])
              if(cor_yh_weight/nfimass>=CorThreshold){
                spe_ref_XCorr_yh_decoy[ncor_yh_decoy]<-cor_yh_weight
                nozero_yh_str_fi_intensitynew3<-log2(sum(yh_str_fi_intensitynew3))#*nnozero_yh_str_fi_intensitynew3
                spe_yh_intensity_decoy[ncor_yh_decoy]<-nozero_yh_str_fi_intensitynew3
              }else{
                spe_ref_XCorr_yh_decoy[ncor_yh_decoy]<-0
                spe_yh_intensity_decoy[ncor_yh_decoy]<-0
              }
              #*(1-exp(-(cor_yh_noise_sum-ncor_yh_noise)))
              
            }
          }
          ncor_yh_decoy<-ncor_yh_decoy+1
        }
        spe_ref_XCorr_yh_sum_decoy<-sum(spe_ref_XCorr_yh_decoy)
        spe_yh_intensity_sum_decoy<-sum(spe_yh_intensity_decoy)
        
        ####
        spe_yh_num_decoy<-sum(spe_ref_XCorr_yh_decoy!=0)
        
        if(max(spe_ref_XCorr_yh_decoy)>0){
          max_yh_cor_index<-which.max(spe_ref_XCorr_yh_decoy)
          if(max_yh_cor_index<=(nfi_yhmass-nfimass)/3){
            max_cor_yh_decoy<-ylabel1[nby-nphk-max_yh_cor_index+1]
          }
          else if(max_yh_cor_index<=2*(nfi_yhmass-nfimass)/3 & max_yh_cor_index>(nfi_yhmass-nfimass)/3){
            max_yh_cor_index1<-max_yh_cor_index-(nfi_yhmass-nfimass)/3
            max_cor_yh_decoy<-ylabel2[nby-nphk-max_yh_cor_index1+1]
          }
          else{
            max_yh_cor_index2<-max_yh_cor_index-2*(nfi_yhmass-nfimass)/3
            max_cor_yh_decoy<-ylabel3[nby-nphk-max_yh_cor_index2+1]
          }
        }else{
          max_cor_yh_decoy<-"NA"
        }
        
        ####
        #
        spe_ref_XCorr_mass_bhyf<-c(spe_ref_XCorr_bh,spe_ref_XCorr_yf)
        spe_ref_XCorr_mass_bfyh<-c(spe_ref_XCorr_bf,spe_ref_XCorr_yh)
        spe_ref_XCorr_mass_bhyf_sum<-sum(spe_ref_XCorr_mass_bhyf)
        spe_ref_XCorr_mass_bfyh_sum<-sum(spe_ref_XCorr_mass_bfyh)
        spe_ref_XCorr_mass_bhyf_max<-max(spe_ref_XCorr_mass_bhyf)
        spe_ref_XCorr_mass_bfyh_max<-max(spe_ref_XCorr_mass_bfyh)
        spe_ref_XCorr_mass_bhyf_max_index<-which.max(c(spe_ref_XCorr_bh_sum,spe_ref_XCorr_yf_sum))
        spe_ref_XCorr_mass_bfyh_max_index<-which.max(c(spe_ref_XCorr_bf_sum,spe_ref_XCorr_yh_sum))
        
        spe_num_mass_bhyf_sum<-spe_bh_num+spe_yf_num
        spe_num_mass_bfyh_sum<-spe_bf_num+spe_yh_num
        spe_intensity_mass_bhyf_sum<-spe_bh_intensity_sum+spe_yf_intensity_sum
        spe_intensity_mass_bfyh_sum<-spe_bf_intensity_sum+spe_yh_intensity_sum
        
        mass_bhyf_label<-c(max_cor_bh,max_cor_yf)
        mass_bfyh_label<-c(max_cor_bf,max_cor_yh)
        if(spe_ref_XCorr_mass_bhyf_max>0 & spe_ref_XCorr_mass_bfyh_max>0){
          spe_ref_XCorr_sum_ph[phk]<-paste(spe_ref_XCorr_mass_bhyf_sum,spe_ref_XCorr_mass_bfyh_sum,sep="_")
          spe_fi_num_ph[phk]<-paste(spe_num_mass_bhyf_sum,spe_num_mass_bfyh_sum,sep="_")
          spe_fi_intensity_sum_ph[phk]<-paste(spe_intensity_mass_bhyf_sum,spe_intensity_mass_bfyh_sum,sep="_")
          max_cor_by_ph[phk]<-paste(mass_bhyf_label[spe_ref_XCorr_mass_bhyf_max_index],mass_bfyh_label[spe_ref_XCorr_mass_bfyh_max_index],sep="_")
        }else{
          spe_ref_XCorr_sum_ph[phk]<-"0_0"
          spe_fi_num_ph[phk]<-"0_0"
          spe_fi_intensity_sum_ph[phk]<-"0_0"
          max_cor_by_ph[phk]<-"NA_NA"
        }
        ####
        spe_ref_XCorr_mass_bhyf_decoy<-c(spe_ref_XCorr_bh_decoy,spe_ref_XCorr_yf_decoy)
        spe_ref_XCorr_mass_bfyh_decoy<-c(spe_ref_XCorr_bf_decoy,spe_ref_XCorr_yh_decoy)
        spe_ref_XCorr_mass_bhyf_sum_decoy<-sum(spe_ref_XCorr_mass_bhyf_decoy)
        spe_ref_XCorr_mass_bfyh_sum_decoy<-sum(spe_ref_XCorr_mass_bfyh_decoy)
        spe_ref_XCorr_mass_bhyf_max_decoy<-max(spe_ref_XCorr_mass_bhyf_decoy)
        spe_ref_XCorr_mass_bfyh_max_decoy<-max(spe_ref_XCorr_mass_bfyh_decoy)
        spe_ref_XCorr_mass_bhyf_max_index_decoy<-which.max(c(spe_ref_XCorr_bh_sum_decoy,spe_ref_XCorr_yf_sum_decoy))
        spe_ref_XCorr_mass_bfyh_max_index_decoy<-which.max(c(spe_ref_XCorr_bf_sum_decoy,spe_ref_XCorr_yh_sum_decoy))
        
        spe_num_mass_bhyf_sum_decoy<-spe_bh_num+spe_yf_num_decoy
        spe_num_mass_bfyh_sum_decoy<-spe_bf_num+spe_yh_num_decoy
        spe_intensity_mass_bhyf_sum_decoy<-spe_bh_intensity_sum_decoy+spe_yf_intensity_sum_decoy
        spe_intensity_mass_bfyh_sum_decoy<-spe_bf_intensity_sum_decoy+spe_yh_intensity_sum_decoy
        
        mass_bhyf_label_decoy<-c(max_cor_bh_decoy,max_cor_yf_decoy)
        mass_bfyh_label_decoy<-c(max_cor_bf_decoy,max_cor_yh_decoy)
        if(spe_ref_XCorr_mass_bhyf_max_decoy>0 & spe_ref_XCorr_mass_bfyh_max_decoy>0){
          spe_ref_XCorr_sum_ph_decoy[phk]<-paste(spe_ref_XCorr_mass_bhyf_sum_decoy,spe_ref_XCorr_mass_bfyh_sum_decoy,sep="_")
          spe_fi_num_ph_decoy[phk]<-paste(spe_num_mass_bhyf_sum_decoy,spe_num_mass_bfyh_sum_decoy,sep="_")
          spe_fi_intensity_sum_ph_decoy[phk]<-paste(spe_intensity_mass_bhyf_sum_decoy,spe_intensity_mass_bfyh_sum_decoy,sep="_")
          max_cor_by_ph_decoy[phk]<-paste(mass_bhyf_label_decoy[spe_ref_XCorr_mass_bhyf_max_index_decoy],mass_bfyh_label_decoy[spe_ref_XCorr_mass_bfyh_max_index_decoy],sep="_")
        }else{
          spe_ref_XCorr_sum_ph_decoy[phk]<-"0_0"
          spe_fi_num_ph_decoy[phk]<-"0_0"
          spe_fi_intensity_sum_ph_decoy[phk]<-"0_0"
          max_cor_by_ph_decoy[phk]<-"NA_NA"
        }
      }
    }
    
    spe_ref_XCorr_sum[i]<-paste(spe_ref_XCorr_sum_ph,collapse=";")
    spe_fi_num[i]<-paste(spe_fi_num_ph,collapse=";")
    spe_fi_intensity_sum[i]<-paste(spe_fi_intensity_sum_ph,collapse=";")
    max_cor_by[i]<-paste(max_cor_by_ph,collapse=";")
    ####
    spe_ref_XCorr_sum[i+nms]<-paste(spe_ref_XCorr_sum_ph_decoy,collapse=";")
    spe_fi_num[i+nms]<-paste(spe_fi_num_ph_decoy,collapse=";")
    spe_fi_intensity_sum[i+nms]<-paste(spe_fi_intensity_sum_ph_decoy,collapse=";")
    max_cor_by[i+nms]<-paste(max_cor_by_ph_decoy,collapse=";")
    
    setTxtProgressBar(pb,i)
  }
  
  newdata_phsty<-data_phsty
  newdata_phsty$Cosine_by_Label<-max_cor_by
  newdata_phsty$Spe_Ref_XCorr_Sum<-spe_ref_XCorr_sum
  newdata_phsty$Spe_FragmentIons_Number<-spe_fi_num
  newdata_phsty$Spe_FragmentIons_Intensity_Sum<-spe_fi_intensity_sum
  write.csv(newdata_phsty,file=FinalFileName,row.names=F)
  
  close(pb)
}

#####07. Clear up the phosphosite results from EDIA

PHOSITE<-function(AfterEDIAFileName,ResultFileName,FIonsNumber=4){
  Phosphosite_data1<-read.csv(AfterEDIAFileName,head=T,stringsAsFactors=F,sep=",")
  nPhosphosite_data1<-length(Phosphosite_data1[[1]])
  phospho_fions_num_index<-vector()
  ifions_num_index<-1
  for(i in 1:nPhosphosite_data1){
    phospho_fions<-strsplit(Phosphosite_data1$Fragment_Ions[i],";")
    phospho_fions_num<-length(phospho_fions[[1]])
    if(phospho_fions_num>=FIonsNumber){
      phospho_fions_num_index[ifions_num_index]<-i
      ifions_num_index<-ifions_num_index+1
    }
  }
  Phosphosite_data<-Phosphosite_data1[phospho_fions_num_index,]
  nPhosphosite_data<-length(Phosphosite_data[[1]])
  modseq<-Phosphosite_data$Peptide_Modified_Sequence
  modseq<-gsub("\\[\\+42\\]","a",modseq)
  modseq<-gsub("\\[\\+80\\]","p",modseq)
  modseq<-gsub("\\[\\+16\\]","o",modseq)
  modseq<-gsub("\\[\\+57\\]","c",modseq)
  modseq<-gsub("\\[\\+10\\]","r",modseq)
  modseq<-gsub("\\[\\+8\\]","l",modseq)
  modseq_ph<-gsub("a|o|c|r|l","",modseq)
  nomodseq<-gsub("a|o|c|p|r|l","",modseq)
  
  newfname<-vector()
  newpname<-vector()
  newpepseq<-vector()
  newpms<-vector()
  newpositioninpep<-vector()
  newphosphosite_num<-vector()
  newprecursor_mass<-vector()
  newcharge<-vector()
  newfion<-vector()
  newtran_mass<-vector()
  newrt<-vector()
  newstime<-vector()
  newetime<-vector()
  newsne<-vector()
  new_Spe_Ref_XCorr_Total_Sum<-vector()
  new_Spe_Ref_XCorr_Total_Sum_Filter<-vector()
  new_Found_Frigmentions_Sum_Intensity_Log<-vector()
  new_Found_Frigmentions_Sum_Intensity_Log_Filter<-vector()
  new_Number_FragIons_Total<-vector()
  new_Number_FragIons_aboveCorthreshold_Total<-vector()
  new_Found_Frigmentions_Sum_Intensity_Log_Filter_Average<-vector()
  new_Found_Frigmentions_Sum_Intensity_Total_Filter<-vector()
  
  newspe_ref_XCorr_sum<-vector()
  newspe_fi_num<-vector()
  newspe_fi_intensity_sum<-vector()
  newcosine_label<-vector()
  
  iphospho<-1
  pb<- txtProgressBar(min=1,max=nPhosphosite_data,char="=>")
  for(i in 1:nPhosphosite_data){
    char_modseq_ph<-gregexpr("p",modseq_ph[i])
    Xcosine_data<-strsplit(Phosphosite_data$Spe_Ref_XCorr_Sum[i],";")
    spe_ion_num_data<-strsplit(Phosphosite_data$Spe_FragmentIons_Number[i],";")
    spe_ion_intensity_data<-strsplit(Phosphosite_data$Spe_FragmentIons_Intensity_Sum[i],";")
    spe_ion_by_label_data<-strsplit(Phosphosite_data$Cosine_by_Label[i],";")
    nchar_modseq_ph<-sum(char_modseq_ph[[1]]>=1)
    for(j in 1:nchar_modseq_ph){
      nphk<-char_modseq_ph[[1]][j]-j
      newfname[iphospho]<-Phosphosite_data$File_Name[i]
      newpname[iphospho]<-Phosphosite_data$Protein_Name[i]
      newpepseq[iphospho]<-nomodseq[i]
      newpms[iphospho]<-Phosphosite_data$Peptide_Modified_Sequence[i]
      newpositioninpep[iphospho]<-nphk
      newphosphosite_num[iphospho]<-nchar_modseq_ph
      newprecursor_mass[iphospho]<-Phosphosite_data$Precursor_Mass[i]
      newcharge[iphospho]<-Phosphosite_data$Precursor_Charge[i]
      newfion[iphospho]<-Phosphosite_data$Fragment_Ions[i]
      newtran_mass[iphospho]<-Phosphosite_data$FIMasses[i]
      newrt[iphospho]<-Phosphosite_data$Rentention_Time[i]
      newstime[iphospho]<-Phosphosite_data$Start_Time[i]
      newetime[iphospho]<-Phosphosite_data$End_Time[i]
      newsne[iphospho]<-Phosphosite_data$Serial_Number_Extension[i]
      new_Spe_Ref_XCorr_Total_Sum[iphospho]<-Phosphosite_data$Spe_Ref_XCorr_Total_Sum[i]
      new_Spe_Ref_XCorr_Total_Sum_Filter[iphospho]<-Phosphosite_data$Spe_Ref_XCorr_Total_Sum_Filter[i]
      new_Found_Frigmentions_Sum_Intensity_Log[iphospho]<-Phosphosite_data$Found_Frigmentions_Sum_Intensity_Log[i]
      new_Found_Frigmentions_Sum_Intensity_Log_Filter[iphospho]<-Phosphosite_data$Found_Frigmentions_Sum_Intensity_Log_Filter[i]
      new_Number_FragIons_Total[iphospho]<-Phosphosite_data$Number_FragIons_Total[i]
      new_Number_FragIons_aboveCorthreshold_Total[iphospho]<-Phosphosite_data$Number_FragIons_aboveCorthreshold_Total[i]
      new_Found_Frigmentions_Sum_Intensity_Log_Filter_Average[iphospho]<-Phosphosite_data$Found_Frigmentions_Sum_Intensity_Log_Filter_Average[i]
      new_Found_Frigmentions_Sum_Intensity_Total_Filter[iphospho]<-Phosphosite_data$Found_Frigmentions_Sum_Intensity_Total_Filter[i]
      
      newspe_ref_XCorr_sum1<-Xcosine_data[[1]][j]
      newspe_ref_XCorr_sum1_add<-strsplit(newspe_ref_XCorr_sum1[[1]],"_")
      newspe_ref_XCorr_sum1_add1<-as.numeric(newspe_ref_XCorr_sum1_add[[1]])
      newspe_ref_XCorr_sum[iphospho]<-sum(newspe_ref_XCorr_sum1_add1)
      
      newspe_fi_num1<-spe_ion_num_data[[1]][j]
      newspe_fi_num1_add<-strsplit(newspe_fi_num1[[1]],"_")
      newspe_fi_num1_add1<-as.numeric(newspe_fi_num1_add[[1]])
      newspe_fi_num[iphospho]<-sum(newspe_fi_num1_add1)
      
      newspe_fi_intensity_sum1<-spe_ion_intensity_data[[1]][j]
      newspe_fi_intensity_sum1_add<-strsplit(newspe_fi_intensity_sum1[[1]],"_")
      newspe_fi_intensity_sum1_add1<-as.numeric(newspe_fi_intensity_sum1_add[[1]])
      newspe_fi_intensity_sum[iphospho]<-sum(newspe_fi_intensity_sum1_add1)
      
      newcosine_label[iphospho]<-spe_ion_by_label_data[[1]][j]
      iphospho<-iphospho+1
    }
    setTxtProgressBar(pb, i)
  }
  newisdecoy<-c(rep("Targets",(iphospho-1)/2),rep("Decoys",(iphospho-1)/2))
  #TD_Substract<-newspe_ref_XCorr_sum[1:((iphospho-1)/2)]-newspe_ref_XCorr_sum[((iphospho-1)/2+1):(iphospho-1)]
  #spe_ref_XCorr_Sum_TD_Substract<-c(TD_Substract,-TD_Substract)
  newphospho_dataframe<-data.frame(File_Name=newfname,Protein_Name=newpname,
                                   Peptide_Sequence=newpepseq,Peptide_Modified_Sequence=newpms,
                                   PositionInPeptide=newpositioninpep,PhosphoSite_Number=newphosphosite_num,
                                   Precursor_Mass=newprecursor_mass,
                                   Precursor_Charge=newcharge,Fragment_Ions=newfion,
                                   FIMasses=newtran_mass,
                                   Rentention_Time=newrt,Start_Time=newstime,End_Time=newetime,
                                   Serial_Number_Extension=newsne,
                                   Spe_Ref_XCorr_Total_Sum=new_Spe_Ref_XCorr_Total_Sum,
                                   Spe_Ref_XCorr_Total_Sum_Filter=new_Spe_Ref_XCorr_Total_Sum_Filter,
                                   Found_Frigmentions_Sum_Intensity_Log=new_Found_Frigmentions_Sum_Intensity_Log,
                                   Found_Frigmentions_Sum_Intensity_Log_Filter=new_Found_Frigmentions_Sum_Intensity_Log_Filter,
                                   Number_FragIons_Total=new_Number_FragIons_Total,
                                   Number_FragIons_aboveCorthreshold_Total=new_Number_FragIons_aboveCorthreshold_Total,
                                   Found_Frigmentions_Sum_Intensity_Log_Filter_Average=new_Found_Frigmentions_Sum_Intensity_Log_Filter_Average,
                                   Found_Frigmentions_Sum_Intensity_Total_Filter=new_Found_Frigmentions_Sum_Intensity_Total_Filter,
                                   Cosine_by_Label=newcosine_label,
                                   Spe_Ref_XCorr_Sum=newspe_ref_XCorr_sum,
                                   #Spe_Ref_XCorr_Sum_TD_Substract=spe_ref_XCorr_Sum_TD_Substract,
                                   Spe_FragmentIons_Number=newspe_fi_num,
                                   Spe_FragmentIons_Intensity_Sum=newspe_fi_intensity_sum,
                                   Class=newisdecoy)
  resultname<-paste("PhoResult_total_",ResultFileName,sep="")
  write.csv(newphospho_dataframe,file=resultname,row.names=F)
  
  #####split the zero and the nonzero value
  newphospho_dataframe_decoy_index<-which(newphospho_dataframe$Protein_Name=="Decoys")
  nnewphospho_dataframe_decoy_index<-length(newphospho_dataframe_decoy_index)
  newphospho_dataframe_target<-newphospho_dataframe[-newphospho_dataframe_decoy_index,]
  pho_target_zero<-vector()
  pho_target_zero_index<-1
  for(k in 1:nnewphospho_dataframe_decoy_index){
    target_XCorr_asnumber<-as.numeric(newphospho_dataframe_target$Spe_Ref_XCorr_Sum[k])
    if(target_XCorr_asnumber==0){
      pho_target_zero[pho_target_zero_index]<-k
      pho_target_zero_index<-pho_target_zero_index+1
    }
  }
  pho_decoy_zero<-pho_target_zero+nnewphospho_dataframe_decoy_index
  pho_zero<-c(pho_target_zero,pho_decoy_zero)
  newphospho_dataframe_nozero<-newphospho_dataframe[-pho_zero,]
  photarget_name<-paste("PhoResult_",ResultFileName,sep="")
  write.csv(newphospho_dataframe_nozero,file=photarget_name,row.names=F)
  
  close(pb)
}

######08. ML:machine learning method to classify targets and decoys

PHOSITE_ML<-function(PhoResultFileName,method="SVM"){
  datamsdf<-read.csv(PhoResultFileName,head=T,stringsAsFactors=F,sep=",")
  datams<-datamsdf[,24:27]
  if(method=="SVM"){
    C_sigma_value<-vector()
    fdr_value_test<-vector()
    fnr_value_test<-vector()
    fdr_fnr_test<-vector()
    j_span<-seq(5,25,0.5)
    cval_span<-c(1:20)
    jj<-1
    pb<- txtProgressBar(min=1,max=length(cval_span)*length(j_span),char="=>")
    for(CVal in cval_span){
      for(j in j_span){
        newdatams_classfier_test<-ksvm(Class~.,data=datams,kernel = "rbfdot",
                                       kpar=list(sigma=j),C=CVal,cross=5)
        newdatams_predict_test<-predict(newdatams_classfier_test,datams)
        Predict_Class_test<-as.character(newdatams_predict_test)
        class_table<-table(datamsdf$Class,Predict_Class_test)
        fdr_value_test[jj]<-class_table[1,2]/sum(class_table[,2])
        fnr_value_test[jj]<-class_table[2,1]/sum(class_table[,1])
        C_sigma_value[jj]<-paste(CVal,j,sep=";")
        fdr_fnr_test[jj]<-class_table[1,2]/sum(class_table[,2])+class_table[2,1]/sum(class_table[,1])
        jj<-jj+1
        setTxtProgressBar(pb, j+(CVal-1)*length(j_span))
      }
      #setTxtProgressBar(pb, CVal)
    }
    fdr_fnr_min_index<-which.min(fdr_fnr_test)
    C_sigma_split<-strsplit(C_sigma_value[fdr_fnr_min_index],";")
    cval_value<-as.numeric(C_sigma_split[[1]][1])
    sigma_value<-as.numeric(C_sigma_split[[1]][2])
    newdatams_classfier<-ksvm(Class~.,data=datams,kernel = "rbfdot",kpar=list(sigma=sigma_value),
                              C=cval_value,cross=5)
  }
  newdatams_predict<-predict(newdatams_classfier,datams)
  datamsdf$Predict_Class<-as.character(newdatams_predict)
  #newdatams_predict_probable<-predict(newdatams_classfier,datams,type="probabilities")
  #datamsdf$Predict_Probabilities<-newdatams_predict_probable[,2]
  datamsdf_class_index<-which(datamsdf$Class=="Targets")
  datamsdf_target<-datamsdf[datamsdf_class_index,]
  datamsdf_decoy<-datamsdf[-datamsdf_class_index,]
  phoresult_file_name1<-gsub("\\.csv","",PhoResultFileName)
  final_result_name<-paste(phoresult_file_name1,"_ML_finalresult.csv",sep="")
  write.csv(datamsdf_target,file=final_result_name,row.names=F)
  final_result_name<-paste(phoresult_file_name1,"_ML_TD_finalresult.csv",sep="")
  write.csv(datamsdf,file=final_result_name,row.names=F)
  close(pb)
}

######09. Get the Phosphosites from ETRF results
GPHOSITE_ETRF<-function(MSMSFileName="msms.txt",EvidenceFileName="evidence.txt",
                        OutTotalFileName,frommsms="no"){
  data_out_total<-read.csv(OutTotalFileName,head=T,stringsAsFactors=F,sep=",")
  nPhosphosite_data<-length(data_out_total[[1]])
  modseq_out_total<-data_out_total$Peptide_Modified_Sequence
  modseq_out_total<-gsub("\\[\\+42\\]","a",modseq_out_total)
  modseq_out_total<-gsub("\\[\\+80\\]","p",modseq_out_total)
  modseq_out_total<-gsub("\\[\\+16\\]","o",modseq_out_total)
  modseq_out_total<-gsub("\\[\\+10\\]","r",modseq_out_total)
  modseq_out_total<-gsub("\\[\\+8\\]","l",modseq_out_total)
  modseq_out_total<-gsub("\\[\\+57\\]","",modseq_out_total)
  modseq_ph<-gsub("a|o|r|l","",modseq_out_total)
  nomodseq<-gsub("a|o|p|r|l","",modseq_out_total)
  newfname<-vector()
  newpname<-vector()
  newpepseq<-vector()
  newpms<-vector()
  newpositioninpep<-vector()
  newphonumber<-vector()
  newprecursor_mass<-vector()
  newcharge<-vector()
  newfion<-vector()
  newtran_mass<-vector()
  newarea<-vector()
  newrt<-vector()
  newstime<-vector()
  newetime<-vector()
  iphospho<-1
  for(i in 1:nPhosphosite_data){
    char_modseq_ph<-gregexpr("p",modseq_ph[i])
    nchar_modseq_ph<-sum(char_modseq_ph[[1]]>=1)
    for(j in 1:nchar_modseq_ph){
      nphk<-char_modseq_ph[[1]][j]-j
      newfname[iphospho]<-data_out_total$File_Name[i]
      newpname[iphospho]<-data_out_total$Protein_Name[i]
      newpepseq[iphospho]<-nomodseq[i]
      newpms[iphospho]<-modseq_out_total[i]
      newpositioninpep[iphospho]<-nphk
      newphonumber[iphospho]<-nchar_modseq_ph
      newprecursor_mass[iphospho]<-data_out_total$Precursor_Mass[i]
      newcharge[iphospho]<-data_out_total$Precursor_Charge[i]
      newfion[iphospho]<-data_out_total$Fragment_Ions[i]
      newtran_mass[iphospho]<-data_out_total$FIMasses[i]
      newarea[iphospho]<-data_out_total$Area[i]
      newrt[iphospho]<-data_out_total$Rentention_Time[i]
      newstime[iphospho]<-data_out_total$Start_Time[i]
      newetime[iphospho]<-data_out_total$End_Time[i]
      iphospho<-iphospho+1
    }
  }
  newphosphotxt_out_total<-data.frame(File_Name=newfname,Protein_Name=newpname,
                                      Peptide_Sequence=newpepseq,Peptide_Modified_Sequence=newpms,
                                      PositionInPeptide=newpositioninpep,
                                      PhosphoSite_Number=nchar_modseq_ph,
                                      Precursor_Mass=newprecursor_mass,
                                      Precursor_Charge=newcharge,Fragment_Ions=newfion,
                                      FIMasses=newtran_mass,Area=newarea,
                                      Rentention_Time=newrt,Start_Time=newstime,End_Time=newetime)
  ####
  data_msms<-read.table(MSMSFileName,head=T,stringsAsFactors=F,sep="\t")
  data_evi<-read.table(EvidenceFileName,head=T,stringsAsFactors=F,sep="\t")
  #####filter msms.txt
  contami_num_evi<-which(data_evi$Potential.contaminant=="+")
  contamin_pep<-unique(c(data_evi$Modified.sequence[contami_num_evi]))
  ncontamin_pep<-length(contamin_pep)
  contamin_msms<-vector()
  for(ic in 1:ncontamin_pep){
    contamin_msms1<-which(data_msms$Modified.sequence==contamin_pep[ic])
    contamin_msms[ic]<-paste(contamin_msms1,collapse=";")
  }
  contamin_msms_total1<-paste(contamin_msms,collapse=";")
  contamin_msms_total<-strsplit(contamin_msms_total1,";")
  contamin_msms_totalnum<-unique(c(as.numeric(contamin_msms_total[[1]])))
  phnum<-which(data_msms$Phospho..STY.==0)
  reversenum<-which(data_msms$Reverse=="+")
  
  num<-unique(c(phnum,reversenum,na.omit(contamin_msms_totalnum)))
  newdata_msms_two<-data_msms[-num,]
  nnewdata_msms_two<-length(newdata_msms_two$id)
  mult_mod_index<-vector()
  jmult_mod<-1
  #Remove the multi-modified aas
  for(i in 1:nnewdata_msms_two){
    mult_mod<-gregexpr("^_\\(\\w+\\)\\w{1}\\(",newdata_msms_two$Modified.sequence[i])
    nmult_mod<-sum(mult_mod[[1]]>=1)
    if(nmult_mod>0){
      mult_mod_index[jmult_mod]<-i
      jmult_mod<-jmult_mod+1
    }
  }
  if(length(mult_mod_index)==0){
    newdata_msms_one<-newdata_msms_two
  }else{
    newdata_msms_one<-newdata_msms_two[-mult_mod_index,]
  }
  nomodseq1<-as.character(newdata_msms_one$Sequence)
  modseq1<-as.character(newdata_msms_one$Modified.sequence)
  modseq2<-gsub("_","",modseq1)
  modseq2<-gsub("\\(ac\\)(\\w)","\\1a",modseq2)
  modseq2<-gsub("\\(ox\\)","o",modseq2)
  modseq2<-gsub("\\(ph\\)","p",modseq2)
  modseq2<-gsub("\\(ar\\)","r",modseq2)
  modseq2<-gsub("\\(ly\\)","l",modseq2)
  
  modseq1<-gsub("_|\\(|\\)|ac|ox|ar|ly","",modseq1)
  modseq1<-gsub("ph","p",modseq1)
  seqprob<-as.character(newdata_msms_one$Phospho..STY..Probabilities)
  seqprob<-gsub("\\(|\\)","",seqprob)
  seqsdiff<-as.character(newdata_msms_one$Phospho..STY..Score.Diffs)
  seqsdiff<-gsub("\\(|\\)","",seqsdiff)
  nmodseq1<-length(modseq1)
  newnomodseq<-vector()
  newmodseq<-vector()
  newposition<-vector()
  newphonumber<-vector()
  newprob<-vector()
  newsdiff<-vector()
  iposition<-1
  pb<-txtProgressBar(min=1,max=nmodseq1,char="=>")
  for(i in 1:nmodseq1){
    ph_extract<-str_extract_all(modseq1[i],"[S,T,Y]p{0,1}")
    prob_extract<-str_extract_all(seqprob[i],"[S,T,Y]\\d{0,1}\\.{0,1}\\d{0,5}")
    sdiff_extract<-str_extract_all(seqsdiff[i],"[S,T,Y]-{0,1}\\d+\\.{0,1}\\d*")
    char_ph_extract<-gregexpr("p",ph_extract[[1]])
    char_ph_extract_index<-which(char_ph_extract>0)
    char_modseq_ph<-gregexpr("p",modseq1[i])
    nseqph<-sum(char_modseq_ph[[1]]>=1)
    for(j in 1:nseqph){
      nphk<-char_modseq_ph[[1]][j]-j
      newnomodseq[iposition]<-nomodseq1[i]
      newmodseq[iposition]<-modseq2[i]
      newposition[iposition]<-nphk
      newphonumber[iposition]<-nseqph
      prob_extract_num<-str_extract(prob_extract[[1]][char_ph_extract_index[j]],"\\d+.{0,1}\\d*")
      sdiff_extract_num<-str_extract(sdiff_extract[[1]][char_ph_extract_index[j]],"-{0,1}\\d+.{0,1}\\d*")
      newprob[iposition]<-as.numeric(prob_extract_num)
      newsdiff[iposition]<-as.numeric(sdiff_extract_num)
      iposition<-iposition+1
    }
    setTxtProgressBar(pb,i)
  }
  modseq_position_two<-paste(newmodseq,newposition,sep="_")
  no_modseq1<-gsub("p","",newmodseq)
  no_modseq2<-gsub("o","",no_modseq1)
  no_modseq<-gsub("a","",no_modseq2)
  no_modseq<-gsub("r","",no_modseq)
  no_modseq<-gsub("l","",no_modseq)
  onlyseq_position_two<-paste(no_modseq,newposition,sep="_")
  
  newdata_position_two<-data.frame(Sequence=newnomodseq,Modified_Sequence=newmodseq,
                                   PositionInPep=newposition,
                                   PhosphoSite_Number=newphonumber,
                                   Sequence_Position=modseq_position_two,
                                   Only_Sequence_Position=onlyseq_position_two,
                                   Phosphosite_Probability=newprob,
                                   Phosphosite_Scorediff=newsdiff)
  close(pb)
  if(frommsms=="no"){
    out_total_pms<-as.character(newphosphotxt_out_total$Peptide_Modified_Sequence)
    data_position_char<-as.character(newdata_position_two$Modified_Sequence)
    nout_total_pms<-length(out_total_pms)
    out_total_pms_index<-vector()
    iout_total<-1
    for(i in 1:nout_total_pms){
      data_position_char_index1<-which(data_position_char==out_total_pms[i])
      data_position_char_index<-paste(data_position_char_index1,collapse=";")
      out_total_pms_index[iout_total]<-data_position_char_index
      iout_total<-iout_total+1
    }
    newout_total_pms_index<-paste(out_total_pms_index,collapse=";")
    newout_total_pms_index1<-strsplit(newout_total_pms_index,";")
    newout_total_pms_indexnum<-as.numeric(newout_total_pms_index1[[1]])
    newdata_position<-newdata_position_two[newout_total_pms_indexnum,]
    modseq_position<-as.character(newdata_position$Only_Sequence_Position)
  }else{
    newdata_position<-newdata_position_two
    modseq_position<-as.character(newdata_position$Only_Sequence_Position)
  }
  
  #get the largest probability and scorediff
  table_position<-table(modseq_position)
  ntable_position<-length(table_position)
  newnomodseq1<-vector()
  newmodseq1<-vector()
  newposition1<-vector()
  newphonumber1<-vector()
  newseq_position1<-vector()
  newprob1<-vector()
  newsdiff1<-vector()
  pb1<-txtProgressBar(min=1,max=ntable_position,char="=>")
  for(i in 1:ntable_position){
    pms<-dimnames(table_position)[[1]][i]
    pms_index<-which(as.character(newdata_position$Only_Sequence_Position)==pms)
    npms_index<-length(pms_index)
    newnomodseq1[i]<-as.character(newdata_position$Sequence[pms_index[1]])
    newmodseq1[i]<-as.character(newdata_position$Modified_Sequence[pms_index[1]])
    newposition1[i]<-as.numeric(newdata_position$PositionInPep[pms_index[1]])
    newphonumber1[i]<-as.numeric(newdata_position$PhosphoSite_Number[pms_index[1]])
    newseq_position1[i]<-pms
    fion_prob<-vector()
    fion_sdiff<-vector()
    for(j in 1:npms_index){
      fion_prob[j]<-as.numeric(newdata_position$Phosphosite_Probability[pms_index[j]])
      fion_sdiff[j]<-as.numeric(newdata_position$Phosphosite_Scorediff[pms_index[j]])
    }
    fion_probmaxindex<-which.max(fion_prob)
    newprob1[i]<-max(fion_prob)
    newsdiff1[i]<-fion_sdiff[fion_probmaxindex]
    
    setTxtProgressBar(pb1,i)
  }
  seq_position<-paste(newnomodseq1,newposition1,sep="_")
  newdata_position1<-data.frame(Sequence=newnomodseq1,Modified_Squence=newmodseq1,
                                PositionInPep=newposition1,
                                PhosphoSite_Number=newphonumber1,
                                ModifiedSequence_Position=newseq_position1,
                                Sequence_Position=seq_position,
                                Phosphosite_Probability=newprob1,
                                Phosphosite_Scorediff=newsdiff1)
  if(frommsms=="no"){
    phosite_name<-paste("FromOuttotal_Phosphosite_",OutTotalFileName,sep="")
  }else{
    phosite_name<-paste("Frommsms_Phosphosite_",OutTotalFileName,sep="")
  }
  write.csv(newdata_position1,file=phosite_name,row.names=F)
  close(pb1)
}

######10. Get EDIA Phosphosite after machine learning
GPHOSITE_EDIA<-function(NormalFileName,AfterMLFileName){
  data_finalms2_equal<-read.csv(NormalFileName,head=T,stringsAsFactors=F,sep=",")
  nPhosphosite_data<-length(data_finalms2_equal[[1]])
  modseq1<-data_finalms2_equal$Peptide_Modified_Sequence
  modseq1<-gsub("\\[\\+42\\]","a",modseq1)
  modseq1<-gsub("\\[\\+80\\]","p",modseq1)
  modseq1<-gsub("\\[\\+16\\]","o",modseq1)
  modseq1<-gsub("\\[\\+10\\]","r",modseq1)
  modseq1<-gsub("\\[\\+8\\]","l",modseq1)
  modseq1<-gsub("\\[\\+57\\]","",modseq1)
  modseq_equal<-gsub("\\[.+\\]","",modseq1)
  modseq_ph<-gsub("a|o|r|l","",modseq_equal)
  nomodseq<-gsub("a|o|p|r|l","",modseq_equal)
  newfname<-vector()
  newpname<-vector()
  newpepseq<-vector()
  newpms<-vector()
  newpositioninpep<-vector()
  phosnumber<-vector()
  newprecursor_mass<-vector()
  newcharge<-vector()
  newfion<-vector()
  newtran_mass<-vector()
  newrt<-vector()
  newstime<-vector()
  newetime<-vector()
  newsne<-vector()
  new_Spe_Ref_XCorr_Total_Sum<-vector()
  new_Spe_Ref_XCorr_Total_Sum_Filter<-vector()
  new_Found_Frigmentions_Sum_Intensity_Log<-vector()
  new_Found_Frigmentions_Sum_Intensity_Log_Filter<-vector()
  new_Number_FragIons_Total<-vector()
  new_Number_FragIons_aboveCorthreshold_Total<-vector()
  new_Found_Frigmentions_Sum_Intensity_Log_Filter_Average<-vector()
  new_Found_Frigmentions_Sum_Intensity_Total_Filter<-vector()
  iphospho<-1
  for(i in 1:nPhosphosite_data){
    char_modseq_ph<-gregexpr("p",modseq_ph[i])
    nchar_modseq_ph<-sum(char_modseq_ph[[1]]>=1)
    for(j in 1:nchar_modseq_ph){
      nphk<-char_modseq_ph[[1]][j]-j
      newfname[iphospho]<-data_finalms2_equal$File_Name[i]
      newpname[iphospho]<-data_finalms2_equal$Protein_Name[i]
      newpepseq[iphospho]<-nomodseq[i]
      newpms[iphospho]<-data_finalms2_equal$Peptide_Modified_Sequence[i]
      newpositioninpep[iphospho]<-nphk
      phosnumber[iphospho]<-nchar_modseq_ph
      newprecursor_mass[iphospho]<-data_finalms2_equal$Precursor_Mass[i]
      newcharge[iphospho]<-data_finalms2_equal$Precursor_Charge[i]
      newfion[iphospho]<-data_finalms2_equal$Fragment_Ions[i]
      newtran_mass[iphospho]<-data_finalms2_equal$FIMasses[i]
      newrt[iphospho]<-data_finalms2_equal$Rentention_Time[i]
      newstime[iphospho]<-data_finalms2_equal$Start_Time[i]
      newetime[iphospho]<-data_finalms2_equal$End_Time[i]
      newsne[iphospho]<-data_finalms2_equal$Serial_Number_Extension[i]
      new_Spe_Ref_XCorr_Total_Sum[iphospho]<-data_finalms2_equal$Spe_Ref_XCorr_Total_Sum[i]
      new_Spe_Ref_XCorr_Total_Sum_Filter[iphospho]<-data_finalms2_equal$Spe_Ref_XCorr_Total_Sum_Filter[i]
      new_Found_Frigmentions_Sum_Intensity_Log[iphospho]<-data_finalms2_equal$Found_Frigmentions_Sum_Intensity_Log[i]
      new_Found_Frigmentions_Sum_Intensity_Log_Filter[iphospho]<-data_finalms2_equal$Found_Frigmentions_Sum_Intensity_Log_Filter[i]
      new_Number_FragIons_Total[iphospho]<-data_finalms2_equal$Number_FragIons_Total[i]
      new_Number_FragIons_aboveCorthreshold_Total[iphospho]<-data_finalms2_equal$Number_FragIons_aboveCorthreshold[i]
      new_Found_Frigmentions_Sum_Intensity_Log_Filter_Average[iphospho]<-data_finalms2_equal$Found_Frigmentions_Sum_Intensity_Log_Filter_Average[i]
      new_Found_Frigmentions_Sum_Intensity_Total_Filter[iphospho]<-data_finalms2_equal$Found_Frigmentions_Sum_Intensity_Total_Filter[i]
      iphospho<-iphospho+1
    }
  }
  newphosphotxt_equal<-data.frame(File_Name=newfname,Protein_Name=newpname,
                                  Peptide_Sequence=newpepseq,Peptide_Modified_Sequence=newpms,
                                  PositionInPeptide=newpositioninpep,
                                  PhosphoSite_Number=phosnumber,
                                  Precursor_Mass=newprecursor_mass,
                                  Precursor_Charge=newcharge,Fragment_Ions=newfion,
                                  FIMasses=newtran_mass,
                                  Rentention_Time=newrt,Start_Time=newstime,End_Time=newetime,
                                  Serial_Number_Extension=newsne,
                                  Spe_Ref_XCorr_Total_Sum=new_Spe_Ref_XCorr_Total_Sum,
                                  Spe_Ref_XCorr_Total_Sum_Filter=new_Spe_Ref_XCorr_Total_Sum_Filter,
                                  Found_Frigmentions_Sum_Intensity_Log=new_Found_Frigmentions_Sum_Intensity_Log,
                                  Found_Frigmentions_Sum_Intensity_Log_Filter=new_Found_Frigmentions_Sum_Intensity_Log_Filter,
                                  Number_FragIons_Total=new_Number_FragIons_Total,
                                  Number_FragIons_aboveCorthreshold_Total=new_Number_FragIons_aboveCorthreshold_Total,
                                  Found_Frigmentions_Sum_Intensity_Log_Filter_Average=new_Found_Frigmentions_Sum_Intensity_Log_Filter_Average,
                                  Found_Frigmentions_Sum_Intensity_Total_Filter=new_Found_Frigmentions_Sum_Intensity_Total_Filter)
  ####
  data_finalms2<-read.csv(AfterMLFileName,head=T,stringsAsFactors=F,sep=",")
  predict_target_index<-which(data_finalms2$Predict_Class=="Targets")
  new_data_finalms2<-data_finalms2[predict_target_index,]
  final_modseq<-new_data_finalms2$Peptide_Modified_Sequence
  modseq<-c(final_modseq,as.character(newphosphotxt_equal$Peptide_Modified_Sequence))
  modseq<-gsub("\\[\\+42\\]","a",modseq)
  modseq<-gsub("\\[\\+80\\]","p",modseq)
  modseq<-gsub("\\[\\+16\\]","o",modseq)
  modseq<-gsub("\\[\\+8\\]","l",modseq)
  modseq<-gsub("\\[\\+10\\]","r",modseq)
  modseq<-gsub("\\[\\+57\\]","",modseq)
  newfinal_modseq<-gsub("\\[\\+\\d{2,3}\\]","",modseq)
  final_positioninpep<-data_finalms2$PositionInPeptide[predict_target_index]
  final_phosnumber<-data_finalms2$PhosphoSite_Number[predict_target_index]
  newfinal_positioninpep<-c(final_positioninpep,as.numeric(newphosphotxt_equal$PositionInPeptide))
  newfinal_phosnumber<-c(final_phosnumber,as.numeric(newphosphotxt_equal$PhosphoSite_Number))
  newfinal_nomodseq1<-gsub("a","",newfinal_modseq)
  newfinal_nomodseq2<-gsub("p","",newfinal_nomodseq1)
  newfinal_nomodseq<-gsub("o","",newfinal_nomodseq2)
  newfinal_nomodseq<-gsub("r","",newfinal_nomodseq)
  newfinal_nomodseq<-gsub("l","",newfinal_nomodseq)
  final_nomodpep_position<-paste(newfinal_nomodseq,newfinal_positioninpep,sep="_")
  fi_sum_intensity_total_filter<-c(new_data_finalms2$Found_Frigmentions_Sum_Intensity_Total_Filter,
                                   as.numeric(newphosphotxt_equal$Found_Frigmentions_Sum_Intensity_Total_Filter))
  final_modpep_position_df<-data.frame(Modified_Peptide=newfinal_modseq,
                                       Peptide_position=final_nomodpep_position,
                                       PhosphoSite_Number=newfinal_phosnumber,
                                       Found_Frigmentions_Sum_Intensity_Total_Filter=fi_sum_intensity_total_filter)
  ediaFile_Name<-paste("AfterEDIAML_total_",AfterMLFileName,sep="")
  write.csv(final_modpep_position_df,file=ediaFile_Name,row.names=F)
}