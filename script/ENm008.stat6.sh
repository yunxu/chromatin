#!/bin/bash - 
# Calculate each chain contacts based on the Cut

SampleList=($(echo p1500_d3300_c{840,1140}))
#SampleList=($(echo p600_d3300_c{840,1140}))
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
    echo "EnsemblePath= " $EnsemblePath
    echo "OutPath= " $OutPath

    if [ ! -d $OutPath ]; then install -d $OutPath;  fi
    
    CutList=($(printf "%s\n" ${CutList[@]}| sort -n -r))
    for Cut in ${CutList[@]}; do
        CutDir=$OutPath/c$Cut
        if [[ ! -d $CutDir ]] ; then install -d $CutDir; fi
    done
    MaxLogWeight=$(find $EnsemblePath -name \*.pts -print0 -exec grep "# LogWeight" {} \;| cut -d" " -f 3 | sort -n | tail -n1)

    #MaxLogWeight=326.061
    echo "MaxLogWeight=" $MaxLogWeight
    lsfiles=($(find $EnsemblePath -name *_0000.pts | xargs -n1 basename | awk -F_ '{print $1}' | sort -u ))

    for Prefix in ${lsfiles[@]}; do
      echo $Prefix
      ItemID=$Prefix
      JOBNAME=J$ItemID.$SampleSize.$Cut.q
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

EnsemblePath=$EnsemblePath/$Prefix
OutPath=$OutPath/$Prefix
Prefix=$Prefix
CutList=(${CutList[@]})
MaxLogWeight=$MaxLogWeight
Prefix=$Prefix

R --vanilla --slave <<-RSCRIPT
    EnsemblePath=\"$EnsemblePath/$Prefix\"
    OutPath=\"$OutPath\"
    Prefix=\"$Prefix\"
    CutStr = \"${CutList[@]}\"
    MaxLogWeight = $MaxLogWeight

    Cut.Vec= sort(as.vector(unlist(read.table(textConnection(CutStr), stringsAsFactors=FALSE))))
    
    DstFiles.FullName = dir(EnsemblePath,pattern='*.dst',recursive=T,full.names=T)
    DstFiles = basename(DstFiles.FullName)
#DstFiles = DstFiles[1:5]
    # Sample file
    SampleFile = DstFiles.FullName[1]
    SampleDist = read.table(SampleFile)
    NrowIJ = nrow(SampleDist)
    ContactIJ = SampleDist[,1:2]
    NumFiles = length(DstFiles)

    LogWeight = rep(0,NumFiles)
    DistMatrix = matrix(0,nrow=NrowIJ,ncol=NumFiles)
    for ( i in 1:NumFiles){
        DstFile = DstFiles.FullName[i]
        LogWeight[i] = read.table(DstFile,comment.char='',nrows=1)[,3]
        DistMatrix[,i] = read.table(DstFile)[,3]
    }
    ExpDiffLogWeight = sprintf('%.3e',exp(LogWeight-MaxLogWeight))

    for (CutInd in 1:length(Cut.Vec)){
        InCutMatrix = matrix(0,nrow=NrowIJ,ncol=NumFiles)
        for ( i in 1:NumFiles){
            MatchInd = which(DistMatrix[,i]<Cut.Vec[CutInd])
            InCutMatrix[MatchInd,i] = ExpDiffLogWeight[i]
        }
     
        DestDir = paste(OutPath,'/c',Cut.Vec[CutInd],'/',Prefix,sep='') 
        if (!file_test('-d',DestDir)){dir.create(DestDir)}

        ConFiles = paste(DestDir,'/',sub('.dst','.con',DstFiles),sep='')
        print(ConFiles )
        for ( i in 1:NumFiles){
            Comment = paste('# ExpDiffLogWeight= ',ExpDiffLogWeight[i],sep='')
            write.table(file=ConFiles[i],Comment, quote=FALSE,row.names=FALSE,col.names=FALSE)
            write.table(file=ConFiles[i],cbind(ContactIJ,InCutMatrix[,i]), quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
        }
    }
RSCRIPT


    " > $JOBNAME
       qsub $JOBNAME
       rm -f $JOBNAME
      done
done
