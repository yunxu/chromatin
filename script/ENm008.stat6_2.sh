#!/bin/bash - 
# Calculate each chain contacts based on the Cut

SampleList=($(echo p600_d3300_c{840,1140}))
SampleList=($(echo p1500_d3300_c{840,1140}))
SampleSizeList=($(printf "%s\n" ${SampleList[@]} | cut -d_ -f1-2 | sort -u))

for SampleSize in ${SampleSizeList[@]}; do
    echo $SampleSize;
    CutList=()
    for Sample in ${SampleList[@]};do
        if [[ $Sample =~ $SampleSize ]]; then
            Cut=$(echo $Sample | awk -F"_c" '{print $2}')
            CutList=(${CutList[@]} $Cut)
        fi
    done

    EnsemblePath=../result/ENm008/ensemble/$SampleSize
    OutPath=../result/analysis/ENm008/cut/$SampleSize 

    if [ ! -d $OutPath ]; then install -d $OutPath;  fi
    
    CutList=($(printf "%s\n" ${CutList[@]}| sort -n -r))
    for Cut in ${CutList[@]}; do
        CutDir=$OutPath/c$Cut
        if [[ ! -d $CutDir ]] ; then install -d $CutDir; fi
    done
    OutFile=tmp/$SampleSize.comb.c.txt
      ItemID=$SampleSize
      JOBNAME=J.$ItemID.$SampleSize.$Cut.q
      FileLog=$0.log
      FJOB=$JOBNAME

      echo "
#!/bin/bash -
#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y
#$ -M yunxu
#$ -m n
#$ -cwd

EnsemblePath=$EnsemblePath
OutPath=$OutPath
Prefix=$Prefix
CutList=(${CutList[@]})
Prefix=$Prefix

R --vanilla --slave <<-RSCRIPT
    EnsemblePath=\"$EnsemblePath\"
    OutPath=\"$OutPath\"
    Prefix=\"$Prefix\"
    CutStr = \"${CutList[@]}\"
    MaxLogWeight = $MaxLogWeight
    OutFile = '$OutFile'


    Cut.Vec= sort(as.vector(unlist(read.table(textConnection(CutStr), stringsAsFactors=FALSE))))
    
    DstFiles.FullName = dir(EnsemblePath,pattern='*.dst',recursive=T,full.names=T)
# DstFiles.FullName = DstFiles.FullName[1:5]
    # Sample file
    SampleFile = DstFiles.FullName[1]
    SampleDist = read.table(SampleFile)
    NrowIJ = nrow(SampleDist)
    ContactIJ = SampleDist[,1:2]
    NumFiles = length(DstFiles.FullName)

    LogWeight = rep(0,NumFiles)
    DistMatrix = matrix(0,nrow=NrowIJ,ncol=NumFiles)
    for ( i in 1:NumFiles){
        DstFile = DstFiles.FullName[i]
        LogWeight[i] = read.table(DstFile,comment.char='',nrows=1)[,3]
        DistMatrix[,i] = read.table(DstFile)[,3]
    }
    MaxLogWeight = max(LogWeight)
    ExpDiffLogWeight = exp(LogWeight-MaxLogWeight)

    Con.Comb = matrix(0,nrow=NrowIJ,ncol=length(Cut.Vec))
    for (CutInd in 1:length(Cut.Vec)){
        InCutMatrix = matrix(0,nrow=NrowIJ,ncol=NumFiles)
        for ( i in 1:NumFiles){
            MatchInd = which(DistMatrix[,i]<Cut.Vec[CutInd])
            InCutMatrix[MatchInd,i] = ExpDiffLogWeight[i]
        }
        Con.Comb[,CutInd] = sprintf('%.3e',apply(InCutMatrix,1,sum) )  
    }
    CommentStr = paste('# i j',paste(Cut.Vec,collapse=' '), sep=' ')                           
    write.table(file=OutFile,CommentStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(file=OutFile,cbind(ContactIJ,Con.Comb),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep='\t')
RSCRIPT

    " > $JOBNAME
       qsub $JOBNAME
       rm -f $JOBNAME
done
