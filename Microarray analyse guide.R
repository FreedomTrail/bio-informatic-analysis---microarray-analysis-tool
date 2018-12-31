#============================================================================================================================================================================#

#先進行檔案安裝與執行
remove(list=ls())
#must put R code inside analyse file
#excel must be same name as image,negative control chip need surfix "-NC"
#Image must export from genepix pro,with bright , condense and color setting when optimized
setwd("pathway")
source("Microarrayworkv5.R")
#source("Microarrayworkv5.R",encoding="utf-8")
instalmypackage()

#============================================================================================================================================================================#


#analyse data seprately (need run first)
#data type
thresholdmethod<-"positive sample"
#thresholdmethod<-"blank"
proteinconclight<-"650"
avgcutoff<-TRUE
avgcutoffvalue<-100
cvcutoff<-FALSE
cvcutoffvalue<-0.5
ratiocutoff<-TRUE
ratiocutoffvalue<-0.5
thresholdvalue<-c(0,1,2,3)

#決定閥值的判定方法為以sample平均 加上標準差，或是以Blank 當negative control,以LOD當閥值計算，如果是前者：thresholdmethod<-"positive sample"如果是後者，thresholdmethod<-"blank"

Analysedataforseperateresult(thresholdmethod,thresholdvalue,proteinconclight,avgcutoff,avgcutoffvalue,cvcutoff,cvcutoffvalue,ratiocutoff,ratiocutoffvalue)


#============================================================================================================================================================================#

#將data normalized來分析


filetype<-"Mutiple"
proteinconclight<-650
normalizedbasedon<-"avg"
thresholdvalue<-c(0,1,2,3)
thresholdmethod<-"Postive"
preNCavgcutoff<-FALSE
preavgcutoff<-FALSE
preavgcutoffvalue<-100
avgcutoffvalue<-100

precvcutoff<-TRUE
precvcutoffvalue<-0.5

cvcutoff<-FALSE
cvcutoffvalue<-0.5

preratiocutoff<-FALSE
preratiocutoffvalue<-0.5

ratiocutoff<-TRUE
ratiocutoffvalue<-1.0

chipspotvaluetype<-"Raw"

Analysedatanormalize(filetype,proteinconclight,normalizedbasedon,thresholdvalue,thresholdmethod,avgcutoffvalue,cvcutoff,cvcutoffvalue,ratiocutoff,ratiocutoffvalue,precvcutoff,precvcutoffvalue,preratiocutoff,preratiocutoffvalue,preNCavgcutoff,preavgcutoff,preavgcutoffvalue,chipspotvaluetype)

#============================================================================================================================================================================#

#用來新增、修改照片，並測試照片條件
#fortest factor if for check image output format. If TRUE, mean it only generate one hit and information image.
fortest<-FALSE
#if outputdirection is filled it will only generate file in specific folder,if not,it will put in defalt folder
outputdirection<-""
#analysesheetname the sheetname for generate image file
analysesheetname<-"Fake hit"
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

#產生image information
imageinformation<-c("20160825-S","20160901-S","20160902-S","20170221-NC","20170303-01-NC","20170303-02-NC")
showinformation<-FALSE


hitimagecreator(fortest,outputdirection,analysesheetname,hitcountregion,showinformation,imageinformation,intensityranking,hitnamecex,arrangeorder,spotsurroundunglayer,imagearrayrow,imagearraycol,IDsize,IDplace,omaboard,IDline,spotmark,showsampleIDname)

#============================================================================================================================================================================#
#用來合併照片、輸出統整報告
#outputasummaryreport if TRUE
Outputreport<-FALSE
#if want output in specific folder,put in direction
filedirection<-""
#can input data as reference
Inputdatadir<-"2017-05-19 3 chip normalized result.xlsx"
#sheetname for inputdata
Inputdatasheetname<-"Intensitythreshold+2SD"
#deside need output combine image or not
Outputimage<-TRUE

arrangeorder<-"Col"

#row for combine image
imagearrayrow<-3

#col for combine image
imagearraycol<-2
Eyebrowsingresult(Outputreport,filedirection,Inputdatadir,Inputdatasheetname,Outputimage,arrangeorder,imagearrayrow,imagearraycol)

#============================================================================================================================================================================#

       dir.create("new hit", showWarnings = FALSE)

for(filecount in 1:length(diff))
   {
       A<-paste("pathway_1",diff[filecount],sep="")
       B<-paste("pathway_2",diff[filecount],sep="")
       file.copy(A,B)
   }




   outputdirection<-paste(getwd(),"/","Multiple Analyse",sep="")



#============================================================================================================================================================================#


fromdir<-"pathway_1"
todir<-"pathway_2"
   FileImageList<-grep(".j",list.files("pathway_1"),value=T)
   FileImageList2<-grep(".j",list.files("pathway_2"),value=T)
   for( x in 1:length(FileImageList2))
      {
           from<-paste(fromdir,"/",FileImageList2[x],sep="")
           to<-paste(todir,"/",FileImageList2[x],sep="")
           file.copy(from, to, overwrite =TRUE)
      }
   file.copy(from, to, overwrite =
   filereplace<-setdiff(FileImageList,FileImageList2)
#============================================================================================================================================================================#




    tried02<-vector("list",15)
 tried<-vector("list",10)
 t<-1

 for(x in 1:5){  datacheck[[x]]<-list.files(k[x])   }
 for(x in 1:5)
    {
        if(x+1 <=5)
           {
               for(y in (x+1):5)
                   {
                      tried[[t]]<-intersect(datacheck[[x]],datacheck[[y]])
                      t<-t+1
                   }
           }
        
    }

    tried


    #============================================================================================================================================================================#

for(colcount in 1:length(colnames(searchlist)))
                {
                    searchlist[,colcount]<-as.character(searchlist[,colcount])
                }
             colnames(searchlist)<-c("ID","JWID","Gene names","Protein names","Pathway","Gene ontology (GO)")
             searchlist[,"ID"]<-analyseIDresult[[filecount]][,"ID"]
             if(length(which(is.na(searchlist[,"ID"])))!=0)
                {
                     searchlist<-searchlist[-which(is.na(searchlist[,"ID"])),]
                }
             for(placecount in 1:nrow(searchlist))
                {
                    libraryIDlocation<-which(searchlist[placecount,"ID"]==librarydata[,1:2])
                    if(length(libraryIDlocation)!=0)
                       {
                           targetproteinname<-librarydata[libraryIDlocation,1]
                           if(length(targetproteinname)==0)
                              {
                                   targetproteinname<-searchlist[placecount,"ID"]
                                   searchlist[placecount,"JWID"]<-""
                              }else
                                {
                                     searchlist[placecount,"JWID"]<-librarydata[libraryIDlocation,2]
                                }
                       }else
                          {
                               targetproteinname<-searchlist[placecount,"ID"]
                               searchlist[placecount,"JWID"]<-""
                          }

