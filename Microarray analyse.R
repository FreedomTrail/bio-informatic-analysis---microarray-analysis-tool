
#分析資料用的function,將資料夾內的資料按所條件篩選
#for single analyse
#將data normalized來分析
proteinconclight<-650
normalizedbasedon<-"bg"
thresholdvalue<-c(0,1,2,3)
thresholdmethod<-"Postive"
avgcutoff<-TRUE
avgcutoffvalue<-100
cvcutoff<-FALSE
cvcutoffvalue<-0.5
ratiocutoff<-TRUE
ratiocutoffvalue<-1
Analysedataforseperateresult<-function(thresholdmethod,thresholdvalue,proteinconclight,avgcutoff,avgcutoffvalue,cvcutoff,cvcutoffvalue,ratiocutoff,ratiocutoffvalue)
{
     FileList<-grep(".xlsx",list.files(getwd()),value=T)
     FileList<-grep("~",FileList,value=T,invert=T)
     #拿掉normalized的結果
     FileList<-grep("chip normalized result",FileList,value=T,invert=T)
     NCfile<-grep("NC",FileList,value=T)
     Samplefilelist<-grep("NC",FileList,value=T,invert=T)
     FileList<-c(Samplefilelist,NCfile)
     for(filecount in 1:length(FileList))
        {
            
            assayanalysisrsult<-createWorkbook(creator =username())
            datalistname<-c("Chip assay result","Sample value","Blank value")
            datalist<-vector("list",length(datalistname))
            names(datalist)<-datalistname
            #輸入資料
            filepath<-FileList[filecount]
            #讀取檔案，與製作新的表格
            datalist[["Chip assay result"]]<-readWorkbook(filepath, sheet = 1,startRow = 1,colNames = FALSE,rowNames = FALSE)
            datacolnames<-which(datalist[["Chip assay result"]][,1]=="Block")
            assayreport<-datalist[["Chip assay result"]][(datacolnames+1):nrow(datalist[["Chip assay result"]]),]
            colnames(assayreport)<-datalist[["Chip assay result"]][datacolnames,]
            Backgroundchannelcount<-grep("B",colnames(assayreport),value=T)
            Backgroundchannelcount<-grep("50 Median",Backgroundchannelcount,value=T)
            Backgroundchannelcount<-grep("F",Backgroundchannelcount,value=T,invert=T)
            Lightchannelcount<-grep("50 Median - B",colnames(assayreport),value=T)
            if(proteinconclight!="" & length(Lightchannelcount)==2)
               {
                 sampleposition<-which(Lightchannelcount!=grep(proteinconclight,Lightchannelcount,value=T))
                 backgroundposition<-which(Backgroundchannelcount!=grep(proteinconclight,Backgroundchannelcount,value=T))

               }else
                  {
                      sampleposition<-1
                      backgroundposition<-1
                  }
            Lightchannelcount<-union(Lightchannelcount[sampleposition],Lightchannelcount[-sampleposition])
            Backgroundchannelcount<-union(Backgroundchannelcount[backgroundposition],Backgroundchannelcount[-backgroundposition])

            addWorksheet(assayanalysisrsult, "Chip assay result")
            writeData(assayanalysisrsult, "Chip assay result", datalist[["Chip assay result"]], startCol = 1, startRow = 1, xy = NULL,colNames = FALSE, rowNames = FALSE)
            
            
            #去除landmark(row在18,19,20的都要先去除)，並將ID,FM-B移到新表格
            assayreport[,"Row"]<-as.numeric(assayreport[,"Row"])
            assaywithoutlandmark<-subset(assayreport,assayreport[,"Row"]<18)
            assaywithoutlandmark<-assaywithoutlandmark[,c("ID",Lightchannelcount,Backgroundchannelcount,"Block","Row","Column","X","Y")]
            assaywithoutlandmark$Block<-as.numeric(assaywithoutlandmark$Block)
            assaywithoutlandmark$Row<-as.numeric(assaywithoutlandmark$Row)
            assaywithoutlandmark$Column<-as.numeric(assaywithoutlandmark$Column)
            #將原本二重複的表格重複像一到同一row中
            assaydoublet<-select(assaywithoutlandmark[seq(1,nrow(assaywithoutlandmark),by=2),],ID:Y)
            assaydoublet2<-select(assaywithoutlandmark[seq(2,nrow(assaywithoutlandmark),by=2),],ID:Y)
            dotcalresult<-vector("list",length(Lightchannelcount))
            dotcalresultbg<-vector("list",length(Lightchannelcount))
            for( x in 1:length(Lightchannelcount))
               {

                   dotcalresultbg[[x]]<-cbind(t(t(as.numeric(assaydoublet[,x+2]))),t(t(as.numeric(assaydoublet2[,x+2]))),t(t(seq(1:nrow(assaydoublet2)))))
                   dotcalresultbg[[x]][,3]<-apply(dotcalresultbg[[x]][,1:2], 1, mean)
                   dotcalresultbg[[x]]<-t(t(dotcalresultbg[[x]][,3]))
                   colnames(dotcalresultbg[[x]])<-(Backgroundchannelcount[x])
                   dotcalresult[[x]]<-cbind(t(t(as.numeric(assaydoublet[,x+1]))),t(t(as.numeric(assaydoublet2[,x+1]))),t(t(seq(1:nrow(assaydoublet2)))),t(t(seq(1:nrow(assaydoublet2)))))
                   dotcalresult[[x]][,3]<-apply(dotcalresult[[x]][,1:2], 1, mean)
                   dotcalresult[[x]][,4]<-apply(dotcalresult[[x]][,1:2], 1, sd)/dotcalresult[[x]][,3]
                   dotcalresult[[x]]<-as.data.frame(dotcalresult[[x]])
                   dotcalresult[[x]][,4]<-round(dotcalresult[[x]][,4],4)
                   dotcalresultname<-select(assaydoublet,Block:Column)
                   SCN<-t(t(apply(dotcalresultname, 1, spotcoordinatename)))
                   dotcalresult[[x]]<-cbind(assaydoublet$ID,dotcalresult[[x]],dotcalresultbg[[x]],select(assaydoublet,Block:Y),SCN)
                   
                   if(Lightchannelcount[x]=="F550 Median - B550")
                      {
                          colnames(dotcalresult[[x]])<-c("ID","550Spot01","550Spot02","550Average","550C.V.","B550 Median","Block","Row","Column","X","Y","B-R-C")
                      }else
                         {
                              colnames(dotcalresult[[x]])<-c("ID","650Spot01","650Spot02","650Average","650C.V.","B650 Median","Block","Row","Column","X","Y","B-R-C")
                         }

               }
                   
            if(length(Lightchannelcount)==2)
               {
                    if(proteinconclight!="")
                       {
                         assaydoubletratio<-t(t(round(dotcalresult[[1]][,4]/dotcalresult[[2]][,4],4)))
                         colnames(assaydoubletratio)<-"Ratio"
                         assaydoubletbind<-cbind(dotcalresult[[1]][,1:6],dotcalresult[[2]][,2:6],assaydoubletratio,dotcalresult[[1]][,7:12])
                       }else
                          {
                              assaydoubletbind<-cbind(dotcalresult[[1]][,1:6],dotcalresult[[2]][,2:6],dotcalresult[[1]][,7:12])
                          }
               }else
                  {
                      assaydoubletbind<-dotcalresult[[1]]
                  }                               
            #分離blank 跟positive control
            datalist[["Blank value"]]<-rbind(subset(assaydoubletbind,assaydoubletbind[,"ID"] == "Blank"),subset(assaydoubletbind,assaydoubletbind[,"ID"] == "Blnak"),subset(assaydoubletbind,assaydoubletbind[,"ID"] == "Empty"),subset(assaydoubletbind,assaydoubletbind[,"ID"] == "EMPTY"),subset(assaydoubletbind,assaydoubletbind[,"ID"] == "UNKNOWN"))
            datalist[["Blank value"]]<-datalist[["Blank value"]] [order(datalist[["Blank value"]][,"B-R-C"]),]
            datalist[["Sample value"]]<-setdiff(assaydoubletbind,datalist[["Blank value"]])
            #計算閥值(1SD,2SD,3SD)
            if(thresholdmethod=="blank")
               {
                   arrayaverage<-mean(datalist[["Blank value"]][,4])
                   rraysd<-sd(datalist[["Blank value"]][,4])
               }else
                   {
                       arrayaverage<-mean(datalist[["Sample value"]][,4])
                       arraysd<-sd(datalist[["Sample value"]][,4])
                   }
            cutoffresult<-vector("list",length=length(threshold))
            hitcount<-vector("list",length(threshold))
            thresholdname<-thresholdvalue
            hitcountname<-thresholdvalue
            for(thresholdnamecount in 1:length(thresholdvalue))
               {
                   thresholdname[thresholdnamecount]<-paste("Intensitythreshold","+",thresholdvalue[thresholdnamecount],"SD",sep="")
                   hitcountname[thresholdnamecount]<-paste(thresholdvalue[thresholdnamecount],"SDhit",sep="")
               }
            names(cutoffresult)<-thresholdname
            names(hitcount)<-hitcountname
            #以threshold挑選值到三個表格與以去除太低的訊號強度
            
            threshold<-thresholdvalue*arraysd+arrayaverage
            for(x in 1:length(threshold))
               {
                    cutoffresult[[x]]<-subset(datalist[["Sample value"]],datalist[["Sample value"]][,4]>threshold[x])
                    if(avgcutoff==TRUE)
                       {
                           cutoffresult[[x]]<-subset(cutoffresult[[x]],cutoffresult[[x]][,4]>=avgcutoffvalue)
                       }
                    if(cvcutoff==TRUE)
                       {
                           cutoffresult[[x]]<-subset(cutoffresult[[x]],cutoffresult[[x]][,5]<=cvcutoffvalue)
                       }
                    if(proteinconclight!="" & length(Lightchannelcount)==2 & ratiocutoff==TRUE)
                       {
                           cutoffresult[[x]]<-subset(cutoffresult[[x]],cutoffresult[[x]][,"Ratio"]>=ratiocutoffvalue)
                       }     
                    hitcount[[x]]<-nrow(cutoffresult[[x]])
               }
            datalist<-c(datalist,cutoffresult)
            sheetname<-c("Chip assay result","Sample value","Blank value",thresholdname,"Analyse condition")
            #這裡不能用sequence
            for(sheetnamecount in 2:(length(sheetname)-1))
               {
                       addWorksheet(assayanalysisrsult, sheetname[sheetnamecount] )
                       writeDataTable(assayanalysisrsult, sheetname[sheetnamecount], datalist[[sheetnamecount]], startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = FALSE,tableStyle = "TableStyleLight9")
               }
            #記錄condition
            cvcutoffvalue<-ifelse(cvcutoff==TRUE,cvcutoffvalue,"")
            ratiocutoffvalue<-ifelse(ratiocutoff==TRUE,ratiocutoffvalue,"")
            analyse_condition_name<-c("avgcutoffvalue","CV cutoff","thresholdmethod","arrayaverage","arraysd","Ratio cutoff",thresholdname,hitcountname)
            analyse_condition<-data.frame(avgcutoffvalue,cvcutoffvalue,thresholdmethod,arrayaverage,arraysd,ratiocutoffvalue,t(threshold),t(hitcount))
            analyse_condition<-t(analyse_condition)
            colnames(analyse_condition)<-"Value"
            rownames(analyse_condition)<-analyse_condition_name
            analyse_condition<-as.data.frame(analyse_condition)
            addWorksheet(assayanalysisrsult, "Analyse condition" )
            writeDataTable(assayanalysisrsult, "Analyse condition", analyse_condition, startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = TRUE,tableStyle = "TableStyleLight9",withFilter = FALSE)
            #存檔
            saveWorkbook(assayanalysisrsult, filepath, overwrite = TRUE)
            #清除掉所有暫存資料
            rm(assayanalysisrsult,cutoffresult,datalist,hitcount,dotcalresult,assaydoubletbind,assaydoublet,assaydoublet2)
        }
   
}


#============================================================================================================================================================================#

#將data normalized來分析


filetype<-"Mutiple"
proteinconclight<-650
normalizedbasedon<-"avg"
thresholdvalue<-c(0,1,2,3)
thresholdmethod<-"Postive"
preNCavgcutoff<-TRUE
preavgcutoff<-TRUE
preavgcutoffvalue<-100
avgcutoffvalue<-100

precvcutoff<-TRUE

precvcutoffvalue<-0.5
cvcutoff<-TRUE
cvcutoffvalue<-0.5

preratiocutoff<-FALSE
preratiocutoffvalue<-0.5
ratiocutoff<-TRUE
ratiocutoffvalue<-1

chipspotvaluetype<-"Raw"
#Analysedatanormalize(filetype,proteinconclight,normalizedbasedon,thresholdvalue,thresholdmethod,avgcutoffvalue,cvcutoff,cvcutoffvalue,ratiocutoff,ratiocutoffvalue,precvcutoff,precvcutoffvalue,preratiocutoff,preratiocutoffvalue,preNCavgcutoff,preavgcutoff,preavgcutoffvalue,chipspotvaluetype)

Analysedatanormalize<-function(filetype,proteinconclight,normalizedbasedon,thresholdvalue,thresholdmethod,avgcutoffvalue,cvcutoff,cvcutoffvalue,ratiocutoff,ratiocutoffvalue,precvcutoff,precvcutoffvalue,preratiocutoff,preratiocutoffvalue,preNCavgcutoff,preavgcutoff,preavgcutoffvalue,chipspotvaluetype)
   {
       FileList<-grep(".xlsx",list.files(getwd()),value=T)
       FileList<-grep("~",FileList,value=T,invert=T)
       FileList<-grep("chip normalized result",FileList,value=T,invert=T)
       NCfile<-grep("NC",FileList,value=T)
       Samplefilelist<-grep("NC",FileList,value=T,invert=T)
       FileList<-c(Samplefilelist,NCfile)
       analysedatalist<-vector("list",length(FileList))
       analysedatablanklist<-vector("list",length(FileList))
       #讀入sample 與 negative control result
       #製作輸出表格
       normalizeddata<-createWorkbook(creator =username())
       filedatelist<-strsplit(FileList,split="-",fixed=T)
       filedate<-FileList
       sheetnamelist<-FileList
       cvnamelist<-FileList
       blanknamelist<-FileList
       for(filecount in 1:length(FileList))
          {
              filedate[filecount]<-filedatelist[[filecount]][1]
              chiptype<-ifelse(filecount<=length(Samplefilelist),"-S","-NC")
              sheetnamelist[filecount]<-paste(filecount,",",filedate[filecount],chiptype,sep="")
              blanknamelist[filecount]<-paste(filecount,",",filedate[filecount],"-B",sep="")
              cvnamelist[filecount]<-paste(sheetnamelist[filecount]," cv",sep="")
              #需要手動輸入路徑，下方路徑內容要改，如果需要可以更改平均值與CV
              filepath<-FileList[filecount]
              #讀取檔案
              analysedatalist[[filecount]]     <-readWorkbook(filepath, sheet = "Sample value",startRow = 1,colNames = TRUE,rowNames = FALSE)
              analysedatablanklist[[filecount]]<-readWorkbook(filepath, sheet = "Blank value" ,startRow = 1,colNames = TRUE,rowNames = FALSE)
              #輸出原始資料
              addWorksheet(normalizeddata, sheetnamelist[filecount])
              writeDataTable(normalizeddata, sheetnamelist[filecount], analysedatalist[[filecount]], startCol = 1, startRow = 1, xy = NULL,colNames = TRUE, rowNames = FALSE,tableStyle = "TableStyleLight9",withFilter = TRUE)

          } 
       names(analysedatalist)<-sheetnamelist
       names(analysedatablanklist)<-blanknamelist
       inputefilename<-FileList
       #normalize sample result
       if(filetype=="Single")
          {
              chiptype<-ifelse(length(Samplefilelist)>0,"-S","-NC")
              totalchipresult<-datanormalized(chipspotvaluetype,analysedatalist,inputefilename,precvcutoff,preratiocutoff,normalizedbasedon)
              chipspotcrossmedian<-totalchipblankresult[["chipspotcrossmedian"]]
              globalchipmedian<-as.matrix(totalchipblankresult[["globalmedian"]])
              colnames(globalchipmedian)<-"chipmedian"
          }else
             {
                  #normalized sample chip
                  from<-1
                  to<-length(Samplefilelist)
                  mediancount<-vector("numeric",length=(to-from+1))
                  listcount<-from:to
                  for(x in 1:length(mediancount))
                     {
                         mediancount[x]<-median(analysedatalist[[listcount[x]]][,6])
                     }
                  globalnormalizesdbgmedian<-mediancount
                  totalsamplechipresult<-datanormalized(chipspotvaluetype,analysedatalist,from,to,inputefilename,precvcutoff,preratiocutoff,normalizedbasedon)
                  totalchipsampleblankresult<-datanormalized(chipspotvaluetype,analysedatablanklist,from,to,inputefilename,precvcutoff,preratiocutoff,normalizedbasedon)
                  totalchipresult<-totalsamplechipresult
                  chipspotcrossmedian<-totalchipresult[["chipspotcrossmedian"]]
                  globalchipmedian<-totalchipresult[["globalmedian"]]
                  #normalized NC chip
                  from<-length(Samplefilelist)+1
                  to<-length(FileList)
                  mediancount<-vector("numeric",length=(to-from+1))
                  listcount<-from:to
                  for(x in 1:length(mediancount))
                     {
                         mediancount[x]<-median(analysedatalist[[listcount[x]]][,6])
                     }
                  globalnormalizeNCdbgmedian<-mediancount
                  totalNCchipresult<-datanormalized(chipspotvaluetype,analysedatalist,from,to,inputefilename,precvcutoff,preratiocutoff,normalizedbasedon)
                  totalchipNCblankresult<-datanormalized(chipspotvaluetype,analysedatablanklist,from,to,inputefilename,precvcutoff,preratiocutoff,normalizedbasedon)
                  NCchipspotcrossmedian<-totalNCchipresult[["chipspotcrossmedian"]]
                  globalNCchipmedian<-totalNCchipresult[["globalmedian"]]

                  totalchipresult[[1]]<-cbind(totalchipresult[[1]],totalNCchipresult[[1]][,2:ncol(totalchipresult[[1]])])
                  chipspotcrossmedian<-c(chipspotcrossmedian,NCchipspotcrossmedian)
                  globalchipmedian<-t(as.matrix(c(globalchipmedian,globalNCchipmedian)))
                  colnames(globalchipmedian)<-c("Samplechipmedian","NCchipmedian")
             }
       threshold<-thresholdvalue*totalchipresult[["globalnormalizedsd"]]+ totalchipresult[["globalnormalizedavg"]]
       if(thresholdmethod=="Blank")
          {
              totalchipblankresult<-datanormalized(chipspotvaluetype,analysedatablanklist,from,to,inputefilename,precvcutoff,preratiocutoff,normalizedbasedon)
              avgcutoffvalue<-totalchipblankresult[["globalnormalizedavg"]]
              
          }
       #計算sample/NC的ratio
       if(ratiocutoff==TRUE & filetype!="Single")
          {
              compareratio<-matrix(0,ncol=length(globalnormalizesdbgmedian)*length(globalnormalizeNCdbgmedian),nrow=nrow(totalsamplechipresult[["totalchipresult"]]))
              #以下新增

              negativeNCspotvalueposition<-which(totalNCchipresult[["totalchipspotvalue"]]<0,arr.ind=TRUE)
                     totalNCchipresult2<-totalNCchipresult[["totalchipspotvalue"]]
                     for(nvcount in 1:nrow(negativeNCspotvalueposition))
                        {
                            totalNCchipresult2[negativeNCspotvalueposition[nvcount,1],negativeNCspotvalueposition[nvcount,2]]<-0
                        }
              #至此
              for(y in 0:(length(globalnormalizesdbgmedian)-1))
                 {
                     
                     for(x in 1:length(globalnormalizeNCdbgmedian))
                        {
                            if(normalizedbasedon=="avg")
                               {
                                   #compareratio[,x+(length(globalnormalizeNCdbgmedian))*y]<-t(t(  ((totalsamplechipresult[["totalchipspotvalue"]][,y+1]/globalnormalizesdbgmedian[y+1])-(totalNCchipresult[["totalchipspotvalue"]][,x]/globalnormalizeNCdbgmedian[x]))/(abs(totalNCchipresult[["totalchipspotvalue"]][,x]-0)/abs(globalnormalizeNCdbgmedian[x]))  ))
                                   compareratio[,x+(length(globalnormalizeNCdbgmedian))*y]<-t(t(  (totalsamplechipresult[["totalchipspotvalue"]][,y+1]/globalnormalizesdbgmedian[y+1])/(totalNCchipresult2[,x]/globalnormalizeNCdbgmedian[x])  ))

                               }else
                                  {
                                       #compareratio[,x+(length(globalnormalizeNCdbgmedian))*y]<-t(t(  ((totalsamplechipresult[["totalchipspotvalue"]][,y+1]/totalsamplechipresult[["globalmedian"]])-(totalNCchipresult[["totalchipspotvalue"]][,x]/totalNCchipresult[["globalmedian"]]))/(abs(totalNCchipresult[["totalchipspotvalue"]][,x]-0)/abs(totalNCchipresult[["globalmedian"]]))  ))
                                       compareratio[,x+(length(globalnormalizeNCdbgmedian))*y]<-t(t(  (totalsamplechipresult[["totalchipspotvalue"]][,y+1]/totalsamplechipresult[["globalmedian"]])/(totalNCchipresult2[,x]/totalNCchipresult[["globalmedian"]])  ))

                                  }
                        }
                 }
             
              sc<-c(Inf,-Inf)
              scvalue<-c(100,-100)
              for(sccount in 1:2)
                {
                    sclist<-which(compareratio==sc[sccount],arr.ind=TRUE)
                    if(nrow(sclist)!=0)
                       {
                           for(sclistcount in 1:nrow(sclist))
                              {
                                   compareratio[sclist[sclistcount,1],sclist[sclistcount,2]]<-scvalue[sccount]
                              }

                       }
                }
              sc2list<-which(is.na(compareratio),arr.ind=TRUE)
              if(nrow(sc2list)!=0)
                 {
                     for(sclistcount in 1:nrow(sc2list))
                        {
                             compareratio[sclist[sclistcount,1],sc2list[sclistcount,2]]<--200
                        }

                 }



              compareratio<-t(t(apply(compareratio,1,min)))
              #compareratio<-t(t((totalsamplechipresult[["totalchipresult"]][,"spotaverage"]/totalsamplechipresult[["globalnormalizedavg"]])/(totalNCchipresult[["totalchipresult"]][,"spotaverage"]/totalNCchipresult[["globalnormalizedavg"]])))
              colnames(compareratio)<-"Ratio"
              totalchipresult[["totalchipresult"]]<-cbind(totalchipresult[["totalchipresult"]],compareratio)
          }
       totalchipresult[["totalchipresult"]]<-cbind(totalchipresult[["totalchipresult"]],analysedatalist[[1]][,c("Block","Row","Column","X","Y","B-R-C")])    
       hitcount<-vector("list",length(threshold))
       #以threshold挑選值到三個表格與以去除太低的訊號強度
       cutoffresult<-vector("list",length=length(threshold))
       thresholdname<-thresholdvalue
       hitcountname<-thresholdvalue
       for(thresholdnamecount in 1:length(thresholdvalue))
          {
              thresholdname[thresholdnamecount]<-paste("Intensitythreshold","+",thresholdvalue[thresholdnamecount],"SD",sep="")
              hitcountname[thresholdnamecount]<-paste(thresholdvalue[thresholdnamecount],"SDhit",sep="")
          }
       names(cutoffresult)<-thresholdname
       names(hitcount)<-hitcountname
       for(x in 1:length(threshold))
          {
               cutoffresult[[x]]<-subset(totalchipresult[["totalchipresult"]],totalchipresult[["totalchipresult"]][,"spotaverage"]>=threshold[x] & totalchipresult[["totalchipresult"]][,"spotaverage"]>=avgcutoffvalue)
               
               if(cvcutoff==TRUE)
                  {
                      cutoffresult[[x]]<-subset(cutoffresult[[x]],as.numeric(cutoffresult[[x]][,"spotcv"])<=cvcutoffvalue)
                  }
               if(ratiocutoff==TRUE)
                  {
                      cutoffresult[[x]]<-subset(cutoffresult[[x]],as.numeric(cutoffresult[[x]][,"Ratio"])>=ratiocutoffvalue)
                  }
               #改到這裡
               #preavgcutoff<-FALSE
               if(preavgcutoff==TRUE)
                  {
                      avgcolnames<-colnames(totalchipresult[["totalchipspotvalue"]])
                     for(cvcolcount in 1:length(avgcolnames))
                         {
                             cutoffresult[[x]]<-subset(cutoffresult[[x]],cutoffresult[[x]][,avgcolnames[cvcolcount]]>=preavgcutoffvalue)
                         }
                  }
               #preNCavgcutoff<-FALSE
               if(preNCavgcutoff==TRUE &filetype!="Single")
                  {
                      avgcolnames<-colnames(totalNCchipresult[["totalchipspotvalue"]])
                     for(cvcolcount in 1:length(avgcolnames))
                         {
                             cutoffresult[[x]]<-subset(cutoffresult[[x]],cutoffresult[[x]][,avgcolnames[cvcolcount]]<preavgcutoffvalue)
                         }
                  }
                     
               if(precvcutoff==TRUE)
                  {
                      cvcolnames<-colnames(totalchipresult[["totalchipcvvalue"]])
                     for(cvcolcount in 1:length(cvcolnames))
                         {
                             cutoffresult[[x]]<-subset(cutoffresult[[x]],cutoffresult[[x]][,cvcolnames[cvcolcount]]<=precvcutoffvalue)
                         }
                  }
               if(preratiocutoff==TRUE)
                  {
                      ratiocolames<-colnames(totalchipresult[["totalchipratiovalue"]])
                      for(ratiocolcount in 1:ncol(cutoffresult[[x]]))
                         {
                             cutoffresult[[x]]<-subset(cutoffresult[[x]],as.numeric(cutoffresult[[x]][,ratiocolames[ratiocolcount]])>=preratiocutoffvalue)
                         }
                  } 
               Ranking<-t(t(seq(1:nrow(cutoffresult[[x]]))))
               cutoffresult[[x]]<-cutoffresult[[x]][order(cutoffresult[[x]][,"spotaverage"],decreasing=TRUE),]
               colnames(Ranking)<-"Intensityranking"
               cutoffresult[[x]]<-cbind(cutoffresult[[x]],Ranking)      
               hitcount[[x]]<-nrow(cutoffresult[[x]])
               sheetname<-paste("Intensitythreshold","+",thresholdvalue[x],"SD",sep="")
               addWorksheet(normalizeddata, sheetname)
               writeDataTable(normalizeddata, sheetname, cutoffresult[[x]], startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,tableStyle = "TableStyleLight9",withFilter = TRUE)
          }
             
               #記錄condition

               cvcutoffvalue<-ifelse(cvcutoff==TRUE,cvcutoffvalue,"")
               ratiocutoffvalue<-ifelse(ratiocutoff==TRUE,ratiocutoffvalue,"")
               precvcutoffvalue<-ifelse(precvcutoff==TRUE,precvcutoffvalue,"")
               preratiocutoffvalue<-ifelse(preratiocutoff==TRUE,preratiocutoffvalue,"")
               preavgcutoffvalue<-ifelse(preavgcutoff==TRUE,avgcutoffvalue,"")
               preNCavgcutoffvalue<-ifelse(preNCavgcutoff==TRUE,avgcutoffvalue,"")

               analyse_condition_name<-c("filetype","Normalized median sourse","avg cutoff","preavgcutoff","preNCchipavgcutoff","CV cutoff","pre CV cutoff","thresholdmethod","arrayaverage","arraysd","Ratio cutoff","pre Ratio cutoff",thresholdname,hitcountname,names(chipspotcrossmedian),colnames(globalchipmedian))
               analyse_condition<-data.frame(c(filetype,normalizedbasedon,avgcutoffvalue,preavgcutoffvalue,preNCavgcutoff,cvcutoffvalue,precvcutoffvalue,thresholdmethod,totalchipresult[["globalnormalizedavg"]],totalchipresult[["globalnormalizedsd"]],ratiocutoffvalue,preratiocutoffvalue,t(threshold),t(hitcount),chipspotcrossmedian,t(globalchipmedian)))
               analyse_condition<-t(analyse_condition)
               rownames(analyse_condition)<-analyse_condition_name
               analyse_condition<-as.data.frame(analyse_condition)
               addWorksheet(normalizeddata, "Analyse condition")
               writeDataTable(normalizeddata, "Analyse condition", analyse_condition, startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = TRUE,tableStyle = "TableStyleLight9",withFilter = FALSE)
               #存檔
               outputname<-paste(Sys.Date()," ",length(Samplefilelist)," chip normalized result.xlsx",sep="")
               saveWorkbook(normalizeddata, outputname, overwrite = TRUE)

   }
#============================================================================================================================================================================#

#用來新增、修改照片，並測試照片條件
#fortest factor if for check image output format. If TRUE, mean it only generate one hit and information image.
fortest<-FALSE
#if outputdirection is filled it will only generate file in specific folder,if not,it will put in defalt folder
outputdirection<-""
#analysesheetname the sheetname for generate image file
analysesheetname<-"Intensitythreshold+2SD"
#intensityranking for put number ranking followed by signal intensity,must put "Intensity ranking" column
intensityranking<-TRUE
#if hitcountregion auto it will do auto hitcount detect,if format is hitcountregion<-c(X,Y) it will generate image from X to Y.
hitcountregion<-"auto"
#hitnamecex if for information image
hitnamecex<-5
IDsize<-6
#IDplace 1,2,3,4 represent show up in buttom leftside,upper,rightside(clockwise) 
IDplace<-3
omaboard<-2
arrangeorder<-"row"
#nearby spot can not be interger
spotsurroundunglayer<-1
imagearrayrow<-2
imagearraycol<-3
IDline<-1.5
#If true mark spot by yellow rectangle 
spotmark<-TRUE
#TRUE for show protein ID name
showsampleIDname<-TRUE
#hitimagecreator(fortest,outputdirection,analysesheetname,hitcountregion,showinformation,imageinformation,intensityranking,hitnamecex,arrangeorder,spotsurroundunglayer,imagearrayrow,imagearraycol,IDsize,IDplace,omaboard,IDline,spotmark,showsampleIDname)

            
#多樣本比較，必須先製作出比較的表格
hitimagecreator<-function(fortest,outputdirection,analysesheetname,hitcountregion,showinformation,imageinformation,intensityranking,hitnamecex,arrangeorder,spotsurroundunglayer,imagearrayrow,imagearraycol,IDsize,IDplace,omaboard,IDline,spotmark,showsampleIDname)
   {
       #table list
       FileList<-grep(".xlsx",list.files(getwd()),value=T)
       FileList<-grep("~",FileList,value=T,invert=T)
       reportdirection<-grep("result",FileList,value=T)
       FileList<-grep("result",FileList,value=T,invert=T)
       NCfile<-grep("NC",FileList,value=T)
       Samplefilelist<-grep("NC",FileList,value=T,invert=T)
       #要將所有chip sheet sample集合起來輸入資料
       placefilelist<-union(Samplefilelist,NCfile)
       analysedata<-vector("list",length(placefilelist))
       #讀取sample list
       for(filecount in 1:length(placefilelist))
          {
             filepath<-placefilelist[filecount]
             #讀取檔案，與製作新的表格
             analysedata[[filecount]]<-readWorkbook(filepath, sheet = "Sample value",startRow = 1,colNames = TRUE,rowNames = FALSE)
          }
       #image list,讀入jpeg，須注意總容量上限好像是46MB
       FileImageList<-grep(".j",list.files(getwd()),value=T)
       NCImagefile<-grep("NC",FileImageList,value=T)
       SampleImagefilelist<-grep("NC",FileImageList,value=T,invert=T)
       Imagefilenamelist<-union(SampleImagefilelist,NCImagefile)
       Imagefilelist<-vector("list",length(Imagefilenamelist))
       for(filecount in 1:length(Imagefilenamelist))
          {
             filepath<-Imagefilenamelist[filecount]
             Imagefilelist[[filecount]]<-readJPEG(filepath, native = FALSE)
          } 

       #判斷創造的資料夾型態      
       if(outputdirection!="")
          {
               outputfolderdirection<-outputdirection
          }else
             {
                 outputfolderdirection<-"Multiple Analyse"
                 dir.create(outputfolderdirection, showWarnings = FALSE)
                 hitkind<-c("True hit","Fake hit")
                 hitkindfolder<-paste(outputfolderdirection,"/",hitkind,sep="")
                 #if(length(reportdirection)==0)
                  #  {
                   #     hitkindfolder<-union(hitkindfolder,paste(hitkindfolder[1],"/",c("Doublet", "Singlet", "Triplet"),sep=""))
                    #}
                 #將以下的資料夾集中創建
                 for(y in 1:length(hitkindfolder))
                    {
                         dir.create(hitkindfolder[y], showWarnings = FALSE)
                    }

             }

       #輸入所有目標的data
       if(length(reportdirection)!=0)
          {
              analyseIDresult<-vector("list",2)
              filepath<-reportdirection[1]
              #讀取檔案，與製作新的表格
              analyseIDresult[[1]]<-readWorkbook(filepath, sheet = analysesheetname,startRow = 1,colNames = TRUE,rowNames = FALSE)
              sheetcolnames<-colnames(analyseIDresult[[1]])
              if(length(which(sheetcolnames=="B-R-C")==1))
                 {
                     searchingbaseoncol<-"B-R-C"
                 }else
                    {
                        searchingbaseoncol<-"ID"
                    }
          }else
             {

                 analyseIDresult<-vector("list",length(Samplefilelist))
                 for(filecount in 1:length(datafilelist))
                    {
                        #讀取sample list (sample value)
                        if(length(Samplefilelist)!=0)
                           {
                             datafilelist<-Samplefilelist
                           }else
                              {
                                  datafilelist<-NCfile
                              }
                        filepath<-datafilelist[filecount]
                        #讀取檔案，與製作新的表格
                        analyseIDresult[[filecount]]<-readWorkbook(filepath, sheet = analysesheetname,startRow = 1,colNames = TRUE,rowNames = FALSE)
                        #合併list與找出最長的list row，這邊可以重寫成merge
                        for(filecount in 2:length(datafilelist))
                           {
                               analyseIDresult[[1]]<-rbind(analyseIDresult[[1]],analyseIDresult[[filecount]])
                           }
                    }
             }
       #創造使用座標軸的BRC搜尋指令
       analyseBRClist<-t(t(unique(analyseIDresult[[1]][,searchingbaseoncol])))
       colnames(analyseBRClist)<-searchingbaseoncol
       #註解function還沒搞懂怎麼用
       #analyseBRClist<-Reduce(function(...) merge(..., all=T,by=searchingbaseoncol), analyseIDresult)
       
       
       #讀取檔案，與製作新的表格
       #輸出說明用圖片
       inneredgesize<-1
       hitcountinone<-length(Imagefilenamelist)

       if(fortest==FALSE)
          {
              lasthitcount<-nrow(analyseBRClist)
          }else
             {
                 lasthitcount<-1
             }
       #for(hitcount in 1:lasthitcount)
       if(length(hitcountregion)==2)
          {
               hitcountstart<-hitcountregion[1]
               hitcountend<-hitcountregion[2]
          }else
             {
                 hitcountstart<-1
                 hitcountend<-lasthitcount
             }
       for(hitcount in hitcountstart:hitcountend)
          {
              Imagelist<-vector("list",hitcountinone)
              spotplace<-which(analyseBRClist[hitcount,searchingbaseoncol]==analysedata[[1]][,searchingbaseoncol])
              for(spotplacecount in 1:length(spotplace))
                 {
                     spotplace02<-spotplace[spotplacecount]
                     for(filecount in 1:length(Imagefilenamelist))
                        {
                            X<-analysedata[[filecount]][spotplace02,"X"]
                            Y<-analysedata[[filecount]][spotplace02,"Y"]
                            spotrange<-spotarea(X,Y,spotsurroundunglayer)
                            Imagelist[[filecount]]<-Imagefilelist[[filecount]][spotrange[1]: spotrange[2],spotrange[3]: spotrange[4],1:3]
                        }
                     
                     if(searchingbaseoncol=="ID")
                        {
                             BRCcodename<-"0-0-0"
                        }else
                           {
                               BRCcodename<-analyseBRClist[hitcount,searchingbaseoncol]     
                           }
                     if(length(spotplace)==1)
                          {
                     
                              Tagname<-analysedata[[filecount]][spotplace02,"ID"]  
                          }else
                              {
                                    Tagname<-paste(Tagname,"-",spotplacecount,sep="")
                               }
                     #還沒寫完，
                     if(intensityranking==TRUE)
                        {
                             stringlength<-nchar(nrow(analyseBRClist))
                             rankingforshowup<-as.character(analyseIDresult[[1]][hitcount,"Intensityranking"])
                             rankingforfilename<-str_pad(rankingforshowup,stringlength,"left",'0')
                             
                             Imagefilename<-paste(rankingforfilename,",",Tagname,",",BRCcodename,sep="")
                             Tagname<-paste(rankingforshowup,"-",analysedata[[filecount]][spotplace02,"ID"],sep="")
                        }else
                           {
                               Imagefilename<-paste(Tagname,",",BRCcodename,sep="")
                               Tagname<-analysedata[[filecount]][spotplace02,"ID"]
                           }
                     Imagefilename<-paste(outputfolderdirection,"/",Imagefilename,".jpeg",sep="")
                     #fill plot with image
                     hitimagearray(showsampleIDname,Imagelist,Tagname,Imagefilename,imagearrayrow,imagearraycol,inneredgesize,hitcountinone,arrangeorder,IDsize,IDplace,omaboard,IDline,spotmark)   
                 }
             
          }
       if(showinformation==TRUE)
          {
               Imagefilename<-paste(outputfolderdirection,"/","Place Information",".jpeg",sep="")
               informationarray(showsampleIDname,imageinformation,Imagefilename,imagearrayrow,imagearraycol,inneredgesize,hitcountinone,hitnamecex,Imagelist,IDsize,IDplace,omaboard,IDline)          

          } 

       analyse_condition_name<-c("fortest","outputdirection","analysesheetname","intensityranking","hitcountregion","hitnamecex","IDsize","IDplace","omaboard","arrangeorder","spotsurroundunglayer","imagearrayrow","imagearraycol","IDline","spotmark","showsampleIDname")
       analyse_condition<-data.frame(fortest,outputdirection,analysesheetname,intensityranking,hitcountregion,hitnamecex,IDsize,IDplace,omaboard,arrangeorder,spotsurroundunglayer,imagearrayrow,imagearraycol,IDline,spotmark,showsampleIDname)
       analyse_condition<-t(analyse_condition)
       rownames(analyse_condition)<-analyse_condition_name
       filename<-paste(outputfolderdirection,"/",Sys.Date()," ","hitimagecreator-setting.csv",sep="")  
       write.csv(analyse_condition,filename)   
     }
#============================================================================================================================================================================#
#用來合併照片、輸出統整報告
#outputasummaryreport if TRUE
Outputreport<-TRUE
#if want output in specific folder,put in direction
filedirection<-""
#can input data as reference
Inputdatadir<-"~/Desktop/2017-05-15 3 chip normalized result.BG-0R-N.xlsx"
#sheetname for inputdata
Inputdatasheetname<-"Intensitythreshold+1SD"
#deside need output combine image or not
Outputimage<-TRUE
arrangeorder<-"Row"
#row for combine image
imagearrayrow<-2

#col for combine image
imagearraycol<-3
#Eyebrowsingresult(Outputreport,filedirection,Inputdatadir,Inputdatasheetname,Outputimage,imagearrayrow,imagearraycol)

Eyebrowsingresult<-function(Outputreport,filedirection,Inputdatadir,Inputdatasheetname,Outputimage,arrangeorder,imagearrayrow,imagearraycol)
   {
       direction<-list.dirs("Multiple Analyse")
       if(filedirection!="")  
          {
               hitkindcountend<-1
               needposition<-1
               direction<-filedirection
               #imagefilename<-c(""     ,       "Fake hit"  ,  "True hit"     ,    "Doublet hit", "Singlet hit", "Triplet hit")
               imagefilename<-c("True hit")
          }else
             {
                 
                 direction<-grep("Hitarray",direction,value=T,invert=T)
                 needposition<-c(3,2)
                 hitkindcountend<-length(needposition)
                 imagefilename<-c("Fake hit"  ,  "True hit")
             }
       

      
       
       #這座表格的function，Sys.getenv("這裡是電腦Username")要改
       assayanalysisrsult<-createWorkbook(creator =username())
       addWorksheet(assayanalysisrsult, "Multiple analyse result")
       imagefilename<-c("Total hit",imagefilename)
       summary<-matrix(0,nrow=length(imagefilename),ncol=1)
       colnames(summary)<-"Value"
       rownames(summary)<-imagefilename
       for(hitkindcount in 1:length(needposition))
          {
              FileList<-grep(".jp",list.files(direction[needposition[hitkindcount]]),value=T)
              FileList<-subset(FileList,FileList!="Place Information.jpeg")
              Imagelist<-vector("list",length(FileList))
              Hitarrayfolder<-paste(direction[needposition[hitkindcount]],"/","Hitarray",sep="")
              dir.create(Hitarrayfolder, showWarnings = FALSE)
              if(length(Imagelist) !=0)
                 {
                     if(Outputimage==TRUE)
                        {
                            for(filecount in 1:length(FileList))
                               {
                                   singlefilename<-paste(direction[needposition[hitkindcount]],"/",FileList[filecount],sep="")
                                   Imagelist[[filecount]] <- readJPEG(singlefilename)
                               }
                        }
                     
                     
                     hitlist1<-(gsub(".jpeg",replacement="" ,FileList ))
                     hitlist2<-strsplit(hitlist1,split=",",fixed=T)
                     
                     hitlist3<-hitlist2[[1]]
                     for(listlehgth in 2:length(hitlist1))
                        {
                            hitlist3<-rbind(hitlist3,hitlist2[[listlehgth]])
                        }
                     if(length(hitlist2[[1]])==3)
                        {
                             colnames(hitlist3)<-c("Intensityranking","ID","B-R-C")
                             hitlist3[,"Intensityranking"]<-t(t(as.numeric(hitlist3[,"Intensityranking"])))

                        }else
                           {
                               colnames(hitlist3)<-c("ID","B-R-C")
                               
                           }
                     rownames(hitlist3)<-seq(1:length(hitlist1))
                     hitlist3[,"ID"]<-t(t(as.character(hitlist3[,"ID"])))
                     hitlist3[,"B-R-C"]<-t(t(as.character(hitlist3[,"B-R-C"])))
                     if(Inputdatadir!="")
                        {
                             Inputdata<-readWorkbook(Inputdatadir, Inputdatasheetname ,startRow = 1,colNames = TRUE,rowNames = FALSE)

                             for(x in 1:nrow(hitlist3))
                                {
                                     if(as.character(hitlist3[x,"B-R-C"])!="0-0-0")
                                        {
                                             searchcol<-"B-R-C"
                                        }else
                                           {
                                               searchcol<-"ID"
                                           }
                                           

                                           searchname<-hitlist3[x,searchcol]
                                           searchlocation<-which(as.character(Inputdata[,searchcol])==searchname)
                                           inputcolname<-colnames(Inputdata)
                                           if(length(searchlocation)!=0)
                                              {
                                                  tempresult<-as.data.frame(matrix("",nrow=length(searchlocation),ncol=ncol(Inputdata)),stringsAsFactors=FALSE)
                                                  colnames(tempresult)<-inputcolname
                                                  for(y in 1:length(searchlocation))
                                                     {
                                                         tempresult[y,]<-as.character(Inputdata[searchlocation[y],])
                                                     }
                                                 
                                              }else
                                                 {
                                                     tempresult<-as.data.frame(matrix("",nrow=1,ncol=ncol(Inputdata)),stringsAsFactors=FALSE)
                                                     colnames(tempresult)<-inputcolname
                                                     tempresult[1,"ID"]<-hitlist3[x,"ID"]
                                                     tempresult[1,"B-R-C"]<-hitlist3[x,"B-R-C"]
                                                     if(ncol(hitlist3)==3)
                                                        {
                                                            tempresult[1,"Intensityranking"]<-hitlist3[x,"Intensityranking"]
                                                        }
                                                     tempresult[,2]<-"Not found"
                                                 }
                                           if(x==1)
                                              {
                                                tempresult2<-tempresult
                                              }else
                                                 {
                                                     tempresult2<-rbind(tempresult2,tempresult)
                                                 }
                                           
                                }
                             colnames(tempresult2)<-colnames(Inputdata)
                             hitlist3<-tempresult2
                             if(length(which(colnames(Inputdata)=="Intensityranking"))==1)
                                {
                                    hitlist3<-hitlist3[order(as.numeric(hitlist3[,"Intensityranking"]),decreasing=FALSE),]
                                    hitlist3[,"Intensityranking"]<-seq(1:nrow(hitlist3))
                                }

                             rownames(hitlist3)<-seq(1:nrow(hitlist3))
                        }
                     Hitlistname<-paste(Hitarrayfolder,"/",Sys.Date(),"-",imagefilename[needposition[hitkindcount]]," list",".csv",sep="")
                     addWorksheet(assayanalysisrsult, imagefilename[needposition[hitkindcount]])
                     summary[needposition[hitkindcount],1]<-length(FileList)
                     writeDataTable(assayanalysisrsult,imagefilename[needposition[hitkindcount]], hitlist3, startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = FALSE, tableStyle = "TableStyleLight9")
                     write.csv(hitlist3,Hitlistname)

                     if(Outputimage==TRUE)
                        {
                            hitarraycount<-floor(length(FileList)/(imagearraycol*imagearrayrow))
                            hitarrayleft<-length(FileList) %% (imagearraycol*imagearrayrow)
                            operatertimelimit<-ceiling(length(FileList)/(imagearraycol*imagearrayrow))
                            operatertime<-1
                            while(operatertime<=operatertimelimit)
                               {    
                                   #去找可以加入文字的function
                                   graphics.off()
                                   h<-dim(Imagelist[[filecount]])[1]
                                   w<-dim(Imagelist[[filecount]])[2]

                                   #fill plot with image
                                          
                                   usr<-c(0,w,0,h)   
                                   #par(mfrow=c(imagearrayrow,imagearraycol),mar=c(0,0,0,0),oma=c(0,8,0,0),xpd = 0,mgp = c(0, 0, 0),bg=0.9)
                                   if(arrangeorder=="Row")
                                      {
                                           par(mfrow=c(imagearrayrow,imagearraycol),mar=c(1,1,1,1),oma=c(0,0,0,0),xpd = 0,mgp = c(0, 0, 0),bg=0.9)
                                      }else
                                         {
                                             par(mfcol=c(imagearrayrow,imagearraycol),mar=c(1,1,1,1),oma=c(0,0,0,0),xpd = 0,mgp = c(0, 0, 0),bg=0.9)
                                         }
                                   
                                   

                                   if(operatertime!=operatertimelimit | hitarrayleft==0)
                                       {
                                           for(x in (1+(operatertime-1)*(imagearraycol*imagearrayrow)):((imagearraycol*imagearrayrow)+(operatertime-1)*(imagearraycol*imagearrayrow)))
                                              {
                                                  plot(c(0, w), c(0, h), type = "n",axes=F, xlab = "", ylab = "" ,xaxs="i",yaxs="i",xaxt='n',yaxt='n')
                                                  rasterImage(Imagelist[[x]],0,0,w,h)
                            
                                              }
                                       }else
                                           {
                                               for(x in (1+(operatertime-1)*(imagearraycol*imagearrayrow)):(hitarrayleft+(operatertime-1)*(imagearraycol*imagearrayrow)))
                                                  {
                                                      plot(c(0, w), c(0, h), type = "n",axes=F, xlab = "", ylab = "" ,xaxs="i",yaxs="i",xaxt='n',yaxt='n')
                                                      rasterImage(Imagelist[[x]],0,0,w,h)
                                                  }
                                           }
                                   
                                   #add text
                                   #close & save image
                                   Imagefilename<-paste(Hitarrayfolder,"/",Sys.Date(),"-",imagefilename[needposition[hitkindcount]]," array",str_pad(operatertime,2,"left",'0'),".jpeg",sep="")
                                   f<-1
                                   #dev.print(jpeg,filename=Imagefilename,width=f*w+70,height=f*h)
                                   dev.print(jpeg,filename=Imagefilename,width=f*w*imagearraycol,height=f*h*imagearrayrow,bg = "white",quality=100)
                                   dev.off()
                                   if((operatertime==operatertimelimit) & hitarrayleft!=0)
                                      {
                                          lasthitarrayimage<-readJPEG(Imagefilename)
                                          h<-dim(lasthitarrayimage)[1]*((ceiling(hitarrayleft/imagearraycol))/imagearrayrow)
                                          w<-dim(lasthitarrayimage)[2]
                                          modifiedhitarrayimage<-lasthitarrayimage[1:h,1:w,1:3]
                                          writeJPEG(modifiedhitarrayimage,Imagefilename,quality =100, bg = "white")
                                      }
                                   operatertime<-operatertime+1
                               }
                        }
                   
                 }else
                     {
                         hitlist3<-data.frame("No result")
                         Hitlistname<-paste(Hitarrayfolder,"/",Sys.Date(),"-",imagefilename[needposition[hitkindcount]]," list",".csv",sep="")
                         addWorksheet(assayanalysisrsult, imagefilename[needposition[hitkindcount]])
                         writeDataTable(assayanalysisrsult, imagefilename[needposition[hitkindcount]], hitlist3, startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = FALSE,tableStyle = "TableStyleLight9")
                     }
          }
       summary[1,1]<-sum(summary[2:nrow(summary),1])
       summary<-as.data.frame(summary)
       writeDataTable(assayanalysisrsult, "Multiple analyse result", summary, startCol = 1, startRow = 1, xy = NULL,colNames = FALSE,rowNames = TRUE,tableStyle = "TableStyleLight9")
       savefilename<-paste("Multiple Analyse","/",Sys.Date(),"-","Analyse result",".xlsx",sep="")
       if(Outputreport==TRUE)
          {
               saveWorkbook(assayanalysisrsult, savefilename, overwrite = TRUE)
          }
       analyse_condition_name<-c("Outputreport","filedirection","Inputdatadir","Inputdatasheetname","Outputimage","imagearrayrow","imagearraycol")
       analyse_condition<-data.frame(Outputreport,filedirection,Inputdatadir,Inputdatasheetname,Outputimage,imagearrayrow,imagearraycol)
       analyse_condition<-t(analyse_condition)
       colnames(analyse_condition)<-"Value"
       rownames(analyse_condition)<-analyse_condition_name  
       filename<-paste(direction[1],"/",Sys.Date()," ","Eyebrowsingresult-setting.csv",sep="")  
       write.csv(analyse_condition,filename)  
   }
#============================================================================================================================================================================#


Chipimagecombine<-function(arrangeorder,imagearrayrow,imagearraycol)
   {
       #table list
    
       FileImageList<-grep(".j",list.files(getwd()),value=T)
       NCImagefile<-grep("NC",FileImageList,value=T)
       SampleImagefilelist<-grep("NC",FileImageList,value=T,invert=T)
       Imagefilenamelist<-union(NCImagefile,SampleImagefilelist)
       Imagefilelist<-vector("list",length(Imagefilenamelist))
       #讀取image
       for(filecount in 1:length(Imagefilenamelist))
          {
             filepath<-Imagefilenamelist[filecount]
             #讀取檔案，與製作新的表格
             Imagefilelist[[filecount]]<-readJPEG(filepath, native = FALSE)
          }
       
       dir.create("Multiple Analyse", showWarnings = FALSE)
       inneredgesize<-1
       graphics.off()

       omasetting<-c(0,0,0,0)
       par(mar=c(0.5,0.5,0.5,0.5))
       if(arrangeorder=="row")
          {
               par(mfrow=c(imagearrayrow,imagearraycol))
          }else
             {
                 par(mfcol=c(imagearrayrow,imagearraycol))
             }

       par(omi=omasetting)
              for(x in 1:(length(Imagefilenamelist)))
            {    
                 h<-dim(Imagefilelist[[x]])[1]
                 w<-dim(Imagefilelist[[x]])[2] 
                 plot(c(0, w+2*(inneredgesize/2)), c(0, h+2*(inneredgesize/2)), type = "n",axes=F, xlab = "", ylab = "",xaxs="i",yaxs="i")
                 rasterImage(Imagefilelist[[x]],(inneredgesize/2),(inneredgesize/2),w+(inneredgesize/2),h+(inneredgesize/2))
            }     

       f<-(1/4)
       Imagefilename<-paste("Multiple Analyse","/","Chip image combine",".jpeg",sep="")
       dev.print(jpeg,filename=Imagefilename,width=f*imagearraycol*(w+2*(inneredgesize/2))+2,height=f*imagearrayrow*(h+2*(inneredgesize/2))+1,bg = "white",quality=100)

       dev.off()          
   } 
#============================================================================================================================================================================#

Chipimagecompress<-function(reportdirection,imageinformation,arrangeorder,imagearrayrow,imagearraycol)
   {
       FileImageList<-grep(".j",list.files(getwd()),value=T)
       NCImagefile<-grep("NC",FileImageList,value=T)
       SampleImagefilelist<-grep("NC",FileImageList,value=T,invert=T)
       Imagefilenamelist<-union(NCImagefile,SampleImagefilelist)
       Imagefilelist<-vector("list",length(Imagefilenamelist))
       #讀取image
       for(filecount in 1:length(Imagefilenamelist))
          {
             filepath<-Imagefilenamelist[filecount]
             #讀取檔案，與製作新的表格
             Imagefilelist[[filecount]]<-readJPEG(filepath, native = FALSE)
          }
       dir.create("Chip compress", showWarnings = FALSE)
       for(x in 1:(length(Imagefilenamelist)))
          {
              graphics.off()
              omasetting<-c(0,0,0,0)
              par(mar=c(0,0,0,0))
              par(omi=omasetting)
              h<-dim(Imagefilelist[[x]])[1]
              w<-dim(Imagefilelist[[x]])[2] 
              plot(c(0, w), c(0, h), type = "n",axes=F, xlab = "", ylab = "",xaxs="i",yaxs="i")
              rasterImage(Imagefilelist[[x]],0,0,w,h)
              f<-(1/4)
              Imagefilename01<-gsub(".jpg",replacement="" ,Imagefilenamelist[x])
              Imagefilename<-paste("Chip compress","/","C-",Imagefilename01,".tiff",sep="")
              dev.print(tiff,filename=Imagefilename,width=f*w,height=f*h,bg = "white")
              dev.off()   
          }                     
   } 
#============================================================================================================================================================================#
#以下是內部所使用的code，通常不需要用到

#============================================================================================================================================================================#
instalmypackage<-function()
   {
       Pkg <- c("fields","spam","maps","dplyr","base","openxlsx","xts","stringr","jpeg","RCurl","XML")
       inst <- Pkg %in% installed.packages()
       if(length(Pkg[!inst]) > 0) install.packages(Pkg[!inst])
       instpackages <- lapply(Pkg, library, character.only=TRUE)
   }
#============================================================================================================================================================================#
username<-function()
   {
       sysinformation<-Sys.getenv()
          if(length(grep("Apple",names(Sys.getenv(names=TRUE))))==1)
             {
                 systemname<-sysinformation[["USER"]]
             }else
                {
                    systemname<-sysinformation[["USERNAME"]]
                }
       systemname<-as.character(systemname)
       return(systemname)
   }  
#============================================================================================================================================================================#

#將資料整理normalized

datanormalized<-function(chipspotvaluetype,datafornormalized,from,to,inputefilename,precvcutoff,preratiocutoff,normalizedbasedon)
   {
        #這只有在指定先用proteinconccutoff才會使用到
        analysedatalist<-datafornormalized
        filedate<-inputefilename
        sheetnamelist<-inputefilename
        cvnamelist<-inputefilename
        Bnamelist<-inputefilename
        blanknamelist<-inputefilename
        filedatelist<-strsplit(inputefilename,split="-",fixed=T)
        
        cvcount<-0
        for(filecount in from:to)
           {
               cvcount<-cvcount+1
               filedate[filecount]<-filedatelist[[filecount]][1]
               chiptype<-ifelse(grepl("-NC",inputefilename[filecount]),"-NC","-S")
               if(length(which(analysedatalist[[filecount]][,1]=="Blank"))>0)
                  {
                       blanknamelist[filecount]<-paste(cvcount,"-",filedate[filecount],"-B",sep="")
                  }else
                     {
                         sheetnamelist[filecount]<-paste(cvcount,"-",filedate[filecount],chiptype,sep="")
                     }
               cvnamelist[filecount]<-paste(sheetnamelist[filecount]," cv",sep="")
               Bnamelist[filecount]<-paste(sheetnamelist[filecount]," bg",sep="")
               #collect overall spot result
               chipspotvalue<-t(t(analysedatalist[[filecount]][,4]))
               colnames(chipspotvalue)<-sheetnamelist[filecount]
               #collect overall cv result
               chipcvvalue<-t(t(analysedatalist[[filecount]][,5]))
               colnames(chipcvvalue)<-cvnamelist[filecount]
               #collect overall Bg result
               chipBvalue<-t(t(as.matrix(analysedatalist[[filecount]][,6])))
               colnames(chipcvvalue)<-cvnamelist[filecount]
               if(cvcount>=2)
                  {
                      totalchipspotvalue<-cbind(totalchipspotvalue,chipspotvalue)
                      totalchipcvvalue<-cbind(totalchipcvvalue,chipcvvalue)
                      totalchipBvalue<-cbind(totalchipBvalue,chipBvalue)
                  }else
                     {
                         totalchipspotvalue<-chipspotvalue
                         totalchipcvvalue<-chipcvvalue
                         totalchipBvalue<-chipBvalue
                     }

           }

        totalchipspotvalue2<-totalchipspotvalue
        #不確定是否可以將資料以normalized的方式釋出，以前用這種方式，可以用於分析hit  (result is cut by avg result)
        #計算overall median，並normalized
        if(normalizedbasedon=="avg")
           {
               #avg median
               globalmedian<-median(totalchipspotvalue)
               spotcrossmedian<-apply(totalchipspotvalue, 2, median)
           }else
              {
                   #background median
                   globalmedian<-median(totalchipBvalue)
                   spotcrossmedian<-apply(totalchipBvalue, 2, median)
              }
        totalchipspotvalue<-totalchipspotvalue*globalmedian/spotcrossmedian
        #計算normalized後的protein spot avg
        spotcrossavg<-t(t(apply(totalchipspotvalue, 1, mean)))
        colnames(spotcrossavg)<-"spotaverage"
        spotcrosscv<-t(t(apply(totalchipspotvalue, 1, sd)/apply(totalchipspotvalue, 1, mean)))
        #計算normalized後的protein spot cv
        colnames(spotcrosscv)<-"spotcv"
        #製作ID colume
        proteinIDname<-t(t(as.character(analysedatalist[[1]][,1])))
        colnames(proteinIDname)<-"ID"
        #判斷是否要將spot avg資料以原始資料呈現
        if(chipspotvaluetype=="Raw")
           {
               totalchipspotvalue<-totalchipspotvalue2
           }
        #將資料column 合併

        totalchipresult<-cbind(proteinIDname,totalchipspotvalue,spotcrossavg,spotcrosscv)
        if(precvcutoff==TRUE)
           {
               totalchipresult<-cbind(proteinIDname,totalchipspotvalue,spotcrossavg,totalchipcvvalue,spotcrosscv)
           }
        cvcount<-0
        if(preratiocutoff==TRUE)
           {
               cvcount<-cvcount+1
               rationamelist<-inputefilename
               for(filecount in from:to)
                  {
                      rationamelist[filecount]<-paste(sheetnamelist[filecount]," ratio",sep="")
                      #collect overall spot ratio result
                      chipratiovalue<-t(t(analysedatalist[[filecount]][,"Ratio"]))
                      colnames(chipratiovalue)<-rationamelist[filecount]
                      if(filecount>=2)
                         {
                             totalchipratiovalue<-cbind(totalchipratiovalue,chipratiovalue)
                         }else
                            {
                                totalchipratiovalue<-chipratiovalue
                            }
                  }
               totalchipresult<-cbind(totalchipresult,totalchipratiovalue)
           }
       globalnormalizedavg<-mean(spotcrossavg)
       globalnormalizedavgmedian<-median(spotcrossavg)
       globalnormalizedsd<-sd(spotcrossavg)
       totalchipresult<-as.data.frame(totalchipresult)
       for(x in 1: ncol(totalchipresult))
          {
              totalchipresult[,x]<-as.character(totalchipresult[,x])
          }
       for(x in 2:ncol(totalchipresult))
          {
              totalchipresult[,x]<-as.numeric(totalchipresult[,x])
          }
       if(preratiocutoff==TRUE)
           {
               outputresult<-list(totalchipresult=totalchipresult,totalchipcvvalue=totalchipcvvalue,totalchipratiovalue=totalchipratiovalue,spotcrossavg=spotcrossavg,globalnormalizedavg=globalnormalizedavg,spotcrosscv=spotcrosscv,globalnormalizedsd=globalnormalizedsd,globalnormalizedavgmedian=globalnormalizedavgmedian,globalmedian=globalmedian,chipspotvaluetype=chipspotvaluetype,chipspotcrossmedian=spotcrossmedian)
           }else
              {
                  outputresult<-list(totalchipresult=totalchipresult,totalchipspotvalue=totalchipspotvalue,totalchipcvvalue=totalchipcvvalue,spotcrossavg=spotcrossavg,globalnormalizedavg=globalnormalizedavg,spotcrosscv=spotcrosscv,globalnormalizedsd=globalnormalizedsd,globalnormalizedavgmedian=globalnormalizedavgmedian,globalmedian=globalmedian,chipspotvaluetype=chipspotvaluetype,chipspotcrossmedian=spotcrossmedian)
              }
       return(outputresult)
   }
#============================================================================================================================================================================#

#輸出spotimage
hitimagearray<-function(showsampleIDname,Imagelist,Tagname,Imagefilename,imagearrayrow,imagearraycol,inneredgesize,hitcountinone,arrangeorder,IDsize,IDplace,omaboard,IDline,spotmark)
   {
       graphics.off()

       omasetting<-c((1/72),(1/72),(1/72),(1/72))

       omasetting[IDplace]<-omaboard
       par(mar=c(0,0,0,0))
       if(arrangeorder=="row")
          {
               par(mfrow=c(imagearrayrow,imagearraycol))
          }else
             {
                 par(mfcol=c(imagearrayrow,imagearraycol))
             }

       if(showsampleIDname==TRUE)
          {
            
               par(omi=omasetting)
                     
          }else
             {
                  par(omi=c(0,(1/72),(1/72),(1/72)))
             }
       
       for(x in 1:(hitcountinone))
          {    
              h<-dim(Imagelist[[x]])[1]
              w<-dim(Imagelist[[x]])[2] 
              plot(c(0, w+2*(inneredgesize/2)), c(0, h+2*(inneredgesize/2)), type = "n",axes=F, xlab = "", ylab = "",xaxs="i",yaxs="i")
              rasterImage(Imagelist[[x]],(inneredgesize/2),(inneredgesize/2),w+(inneredgesize/2),h+(inneredgesize/2))
              if(spotmark==TRUE)
                 {
                      rect((inneredgesize/2)+w/2-27.5, (inneredgesize/2)+h/2-11.5,(inneredgesize/2)+w/2+27.5, (inneredgesize/2)+h/2+11.5, col = c(NA,0),border = "gold", lwd = 6)
                 }
          }
       f<-4     
       if(showsampleIDname==TRUE)
          {
              if(IDplace==2 | IDplace==4)
                 {
                     mtext(Tagname,IDplace,line=IDline, cex=IDsize, col="black", outer=TRUE)
                     dev.print(jpeg,filename=Imagefilename,width=f*(imagearraycol*(w+2*(inneredgesize/2))+omaboard*18+1),height=f*imagearrayrow*(h+2*(inneredgesize/2))+1,bg = "white",quality=100)
                 }else
                    {
                         mtext(Tagname,IDplace,line=IDline, cex=IDsize*1.86, col="black", outer=TRUE)
                         dev.print(jpeg,filename=Imagefilename,width=f*imagearraycol*(w+2*(inneredgesize/2))+2,height=f*(imagearrayrow*(h+2*(inneredgesize/2))+omaboard*18+1),bg = "white",quality=100)
                    }
          }else
              {
                  dev.print(jpeg,filename=Imagefilename,width=f*imagearraycol*(w+2*(inneredgesize/2))+2,height=f*imagearrayrow*(h+2*(inneredgesize/2))+1,bg = "white",quality=100)
              }
       dev.off()
   }
#============================================================================================================================================================================#

#比較可以分為:全部sample(包含全部NC),須讀取所有file跟focuslist(超過2個同種type以上就需要),一個NC與多sample,一個sample與多NC,且NC一定排最後
informationarray<-function(showsampleIDname,imageinformation,Imagefilename,imagearrayrow,imagearraycol,inneredgesize,hitcountinone,hitnamecex,Imagelist,IDsize,IDplace,omaboard,IDline)
   {
       graphics.off()
       omasetting<-c((5/72),(5/72),(5/72),(5/72))
       omasetting[IDplace]<-omaboard
       par(mar=c(0.5,0.5,0.5,0.5))
       if(arrangeorder=="row")
          {
               par(mfrow=c(imagearrayrow,imagearraycol))
          }else
             {
                 par(mfcol=c(imagearrayrow,imagearraycol))
             }

       if(showsampleIDname==TRUE)
              {
                   par(omi=omasetting)
              }else
                 {
                      par(omi=c(0,(5/72),(5/72),(5/72)))
                 }

       for(x in 1:(hitcountinone))
            {
                 h<-dim(Imagelist[[x]])[1]
                 w<-dim(Imagelist[[x]])[2] 
                 plot(c(0, w+2*(inneredgesize/2)), c(0, h+2*(inneredgesize/2)), type = "n",axes=T, xlab = "", ylab = "",xaxs="i",yaxs="i",xaxt='n',yaxt='n')
                 text((w+2*(inneredgesize/2))/2,(h+2*(inneredgesize/2))*((hitnamecex/7)+4)/8,imageinformation[x],1, cex=hitnamecex, col="black",adj=c(0,0))
            }     
       f<-4
       if(showsampleIDname==TRUE)
          {
              if(IDplace==2 | IDplace==4)
                 {
                     mtext("Protein ID",IDplace,line=IDline, cex=IDsize, col="black", outer=TRUE)
                     dev.print(jpeg,filename=Imagefilename,width=f*(imagearraycol*(w+2*(inneredgesize/2))+omaboard*18+1),height=f*(imagearrayrow*(h+2*(inneredgesize/2))+1),bg = "white",quality=100)
                 }else
                    {
                         mtext("Protein ID",IDplace,line=IDline, cex=IDsize*1.86, col="black", outer=TRUE)
                         dev.print(jpeg,filename=Imagefilename,width=f*imagearraycol*(w+2*(inneredgesize/2))+2,height=f*(imagearrayrow*(h+2*(inneredgesize/2))+omaboard*18+1),bg = "white",quality=100)
                    }
           }else
               {
                   dev.print(jpeg,filename=Imagefilename,width=f*imagearraycol*(w+2*(inneredgesize/2))+2,height=f*imagearrayrow*(h+2*(inneredgesize/2))+1,bg = "white",quality=100)
               }
       dev.off()
   }


#============================================================================================================================================================================#

#如果chip種類不同，視情況有可能要調整以下參數
spotarea<-function(X,Y,spotsurroundunglayer)
   {
       X<-as.numeric(X)/10
       Y<-as.numeric(Y)/10
       spotsurroundunglayer<-as.numeric(spotsurroundunglayer)
       xL<-X-16-50*spotsurroundunglayer
       xR<-X+39+50*spotsurroundunglayer
       yU<-Y+12+23*spotsurroundunglayer
       yD<-Y-11-23*spotsurroundunglayer
       return(c(yD,yU,xL,xR))
   }
#============================================================================================================================================================================#

spotcoordinatename<-function(x)
   {
       x<-str_pad(x,2,"left",'0')
       name<-paste(x[1],x[2],x[3],sep="-")
       return(name)
   }


