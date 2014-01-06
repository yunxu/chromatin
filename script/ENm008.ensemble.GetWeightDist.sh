#!/bin/bash - 
# Generate contact node weighted distance

OutDir=tmp
MaxNum=1000
if [ $(uname) = "Darwin" ];
then
  DataDir=/Users/yunxu/workspace/projects/chromatin/result/ENm008/ensemble/p1500_d3300/ADV
  DataDir=/Users/yunxu/workspace/projects/chromatin/result/ENm008/ensemble/p1500_d3300
else
  DataDir=/dump/yunxu/working/projects.new/chromatin/result/ENm008/ensemble/p1500_d3300/ADV
  DataDir=/dump/yunxu/working/projects.new/chromatin/result/ENm008/ensemble/p1500_d3300
fi

#MaxLogWeight=$(
#  find $DataDir -name \*.pts | \
#    xargs grep Log | sort -k3,3n | tail -n1 | \
#    cut -d" " -f3 )
MaxLogWeight=326.061
echo $MaxLogWeight

partfilename=$1
if [ -n $partfilename  ]; then
  OutDir=$OutDir/$partfilename
  mkdir $OutDir
fi

#R --vanilla --slave <<- RSCRIPT
(
cat <<- RSCRIPT
  MaxNum=$MaxNum
  MaxLogWeight=$MaxLogWeight
  OutDir = "$OutDir"
  DataDir = "$DataDir"
  
#  MaxNum = 1000
#  MaxLogWeight=324.678
#  DataDir =
#  "/Users/yunxu/workspace/projects/chromatin/result/ENm008/ensemble/p1500_d3300/ADV"
#  OutDir = "/tmp"
  NumNodes = 54
  NodeNodeList = t(combn(1:NumNodes,2))
  NodeNodeVec = NodeNodeList[which(NodeNodeList[,2]-NodeNodeList[,1]!=1),]

#  FileListFile = dir(path=DataDir,pattern=".pts",recursive=T,full.names=T) 
  FileListFile = read.table("$partfilename")[,1]

  NewFileNames = NULL
  NumCopyVec = NULL
  for (i in 1:length(FileListFile)){
    FileName = as.character(FileListFile[i])
    LogWeight = read.table(FileName,comment.char="",nrows=1)[,3]
    Weight = exp(LogWeight - MaxLogWeight)
    NumCopies = round(MaxNum * Weight)
    if (NumCopies != 0){
      NumCopyVec = c(NumCopyVec,NumCopies)
      NewFileNames = c(NewFileNames, FileName)
    }
    if( i %% 100 == 0){ 
      cat ("$partfilename", i , "/", length(FileListFile),"\n")
      flush.console() 
    }
  }
#  write.table(file="file.ls.txt",cbind(NewFileNames,NumCopyVec),quote=F,row.names=F,col.names=F) 
#  q()

  iNumFile = length(NewFileNames)
  cat("iNumFile", iNumFile,"\n")
  DistMat = matrix(0, nrow=nrow(NodeNodeVec),ncol=iNumFile)

  for (i in 1:iNumFile){
    cat(i,"/",iNumFile,"\n")
    flush.console()
    FileName = NewFileNames[i]
    PTS = read.table(FileName)[,1:3]
    DistMat[,i] = apply(NodeNodeVec,1,function(x){
      dist( rbind(PTS[x[1],], PTS[x[2],]) )
    })
  }

  for( iRow in 1:nrow(NodeNodeVec)){
    FN = paste(OutDir,"/",NodeNodeVec[iRow,1],"_",NodeNodeVec[iRow,2],".txt",sep="")
    write.table(file=FN,cbind(sprintf("%.3f",DistMat[iRow,]),NumCopyVec),quote=F,row.names=F,col.names=F)
  }
RSCRIPT
) > $partfilename.R

ItemID=$partfilename
JOBNAME=J$ItemID.q
FileLog=$0.log
FJOB=$JOBNAME

(cat <<- QSUB
#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y
#$ -M yunxu
#$ -m n
#$ -cwd

Rscript --vanilla --slave $partfilename.R
QSUB
) > $JOBNAME
qsub $JOBNAME
rm -f $JOBNAME

exit
# ---------------
# merge.sh
# #!/bin/bash 
# OutDir=out
# for nodenodefile in `cat filelist.txt`;
# do
#   echo $nodenodefile
#   for DIR in tmp/*; do
#     cat $DIR/$nodenodefile
#   done >  $OutDir/$nodenodefile
# done
