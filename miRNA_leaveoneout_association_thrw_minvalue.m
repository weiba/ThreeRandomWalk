%%%%%%%%%%loading miRNA similarity data
load ./data/miRNA_similarity_function.txt
M_1 = miRNA_similarity_function;
 
%%%%%%%%%%loading disease similarity data
load ./data/disease_function_similarity_matri_1.txt
D_1 = disease_function_similarity_matri_1;
 
%%%%%%%%%%loading drug similarity data
load ./data/EF_EF_Drug_relationship_matrix.txt
E_1 = EF_EF_Drug_relationship_matrix;

%%%%%%%%%%loading interaction data
load ./data/miRNA_disease_matrix.txt
Y= miRNA_disease_matrix;


%%%%%%%%%%loading  EF and Disease interaction data
load ./data/MicRNA_environment_disease_ass_matrix.txt
Y2 = MicRNA_environment_disease_ass_matrix;

%%%%%%%%%%loading  EF and MicRNA interaction data
load ./data/MicRNA_environment_microRNA_ass_matrix.txt
Y3 = MicRNA_environment_microRNA_ass_matrix;

sum_m_1=1./sum(M_1,2);
sum_D_1=1./sum(D_1,2);
sum_E_1=1./sum(E_1,2);
%%%% initialization
M_size=size(M_1);
D_size=size(D_1);
E_size=size(E_1);

%norminliaztion
xM_1=repmat(sum_m_1,1, M_size(1,1));
M_1=M_1.*xM_1;
xD_1=repmat(sum_D_1,1, D_size(1,1));
D_1=D_1.*xD_1;
xE_1=repmat(sum_E_1,1, E_size(1,1));
E_1=E_1.*xE_1;
 

sum_Y2=sum(sum(Y2));
sum_Y3=sum(sum(Y3));
sum_Y2=1/sum_Y2;
sum_Y3=1/sum_Y3;
Y_2=Y2.*sum_Y2;   %ef-disease
Y_21=Y_2';   %micrna*EF
Y_3=Y3.*sum_Y3;    %ef-mirna
Y_31=Y_3';

sumauc=0;
round=1000;
[row,col]=size(Y);
auc=zeros(col,2);
 TPRAVG_round=zeros(row,1);
 FPRAVG_round=zeros(row,1);
 
for r=1:round
% store evaluation results
r
 TPR=zeros(row,col); % row is 271 ,col is 137
 FPR=zeros(row,col);
 TPRAVG=zeros(row,1);
 FPRAVG=zeros(row,1);
 tpAVG=zeros(row,1);
 fpAVG=zeros(row,1);
 
[Index_PositiveRow,Index_PositiveCol]=find(Y(:,:)==1);
totalassociation=length(Index_PositiveRow);
fold=totalassociation/5;
testset=zeros(fold,1);
trainset=zeros(totalassociation-fold,1);

 p=randperm(totalassociation);
testset=p(1:fold) ;
trainset=p(fold+1:totalassociation) ;
Y1=Y;

 
for i=1:fold
%for i=1:size(Y,2)  % Y is  micrna*disease get column number
    Y1(Index_PositiveRow(testset(i)),Index_PositiveCol(testset(i)))=0;
    
end  
    %%% normilization
    sum_Y1=sum(sum(Y1));
    sum_Y1=1/sum_Y1;
    Y_1=Y1.*sum_Y1;  %micrna-disease
    Y_11=Y_1'; %disease-micrna
     
    prediction=ThrRW_new(M_1,D_1,E_1,Y_1,Y_2,Y_3,Y_31);
    result=prediction.Rmd;  %micrna* disease
   % result1=result';
   % [IX_col,R_col]=sort(result1,'descend');  % 按列排序
   % [IX_row,R_row]=sort(result1,2,'descend');  % 按行排序
   
    tmpmatrix=ones(row,col);
    tmpmatrix1=Y;  
   %  for i=1:size(Index_PositiveRow)-fold
  %      tmpmatrix(Index_PositiveRow(trainset(i)),Index_PositiveCol(trainset(i)))=0;
  %      tmpmatrix1(Index_PositiveRow(trainset(i)),Index_PositiveCol(trainset(i)))=0;  % clean data in trainset 
  %      
   %  end
  
    minvalue=min(result(:));
   if(minvalue>=0)  minvalue=0;
   else
       minvalue=minvalue-0.1;
   end
   result1=result;
     for i=1:size(Index_PositiveRow)-fold
        result1(Index_PositiveRow(trainset(i)),Index_PositiveCol(trainset(i)))=minvalue;% set as minvalue if in trainset . So that they will locate at the end after sorting.  ignore training set when evaluate 
        tmpmatrix1(Index_PositiveRow(trainset(i)),Index_PositiveCol(trainset(i)))=0;  % clean data in trainset 
    end
	
  %   result1=result.*tmpmatrix;
    [IX_col,R_col]=sort(result1,'descend');  % 按列排序
    [IX_row,R_row]=sort(result1,2,'descend');  % 按行排序   micrna*disease
   % EvaluationResults=Evaluationassociation(IX_col,R_col,result1,Y,TPR,FPR,tmpmatrix1);
	  EvaluationResults=Evaluationassociation_minvalue(IX_col,R_col,result1,Y,TPR,FPR,tmpmatrix1,minvalue);
   TPR=EvaluationResults.TPR;
   FPR=EvaluationResults.FPR;

 
 rowsum=sum(tmpmatrix1) ; %returns a row vector containing the sum of each colum
 testdisease=0;
  for i=1:col
      if(rowsum(i)>0)
         testdisease=testdisease+1; 
      end
  end
   for i=1:row
   TPRAVG(i)=sum(TPR(i,:))/testdisease;
   FPRAVG(i)=sum(FPR(i,:))/testdisease;
   end
  TPRAVG_round=TPRAVG_round+TPRAVG;
   FPRAVG_round=FPRAVG_round+FPRAVG;
   
   aucavg=CaculateAUC(FPRAVG,TPRAVG);
  sumauc= sumauc+aucavg;
   
  
   for j=1:col
      auc(j,1)=auc(j,1)+CaculateAUC(FPR(:,j),TPR(:,j));
      if(rowsum(j)>0)
          auc(j,2)=auc(j,2)+1;
      end
   end
   
   
end


sumauc=sumauc/round;
sumauc

TPRAVG_round=TPRAVG_round./round;
FPRAVG_round=FPRAVG_round./round;
aucavg_round=CaculateAUC(FPRAVG_round,TPRAVG_round); 
aucavg_round

 auc_disease=zeros(col,1);
for i=1:col
    if(auc(i,2)>0)
     auc_disease(i)=auc(i,1)/auc(i,2);
    end
end;
auc_disease

