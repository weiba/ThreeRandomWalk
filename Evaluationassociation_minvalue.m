function EvaluationResults = Evaluationassociation_minvalue(IX,R,result,standardMatrix,TPR,FPR,tmpmatrix1,minvalue)  % R,standardMatrix micrna*disease
 
[row,col]=size(R);
%row 
%col
% p=sum(standardMatrix); %求列和
 p=sum(tmpmatrix1); % cleaning the association of standardMatrix in trainset
 p_total=sum(p); %求和
 %p_total


for i=1:col
   if( p(i)>0)
%i=diseaseno;
    for th=1 : row
    tp=0;
    fp=0;
     
    for j=1 : th
        index=int32(R(j,i));
         if(result(index,i)>minvalue)  % if the value is minvalue, it is in the trainset. We should ignore it when evaluate the algorithm 
        if(standardMatrix(index,i)>0)
            tp=tp+1;
        else 
            fp=fp+1;
        end
         end
    end
    if(p(i)>0)
     TPR(th,i)=tp/p(i);
    else
     TPR(th,i)=0;   
    end
    if(row-p(i)>0)
     FPR(th,i)=fp/(row-p(i));
    else
        FPR(th,i)=0;
    end
    
    end
   end % for  

    EvaluationResults.TPR=TPR;
    EvaluationResults.FPR=FPR;
    %EvaluationResults.AVGTPR=TPRAVG;
    %EvaluationResults.AVGFPR=FPRAVG;
end