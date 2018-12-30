function state = ThrRW_new(M_1,D_1,E_1,Y_1,Y_2,Y_3,Y_31)
l1=1;
r1=1;
l2=1;
r2=1;
l3=1;
r3=1;
a=0.9;
Red=Y_2;
Rmd=Y_1;
Rem=Y_3;
parmeter=[l1,r1,l2,r2,l3,r3];
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%leave one percent evaluation
for t=1:max(parmeter)
    nm1=0;
    nd1=0;
    ne1=0;
    nm2=0;
    nd2=0;
    ne2=0;
	nem2=0;
	nem1=0;
	nm3=0;
    nd3=0;
    ne3=0;
	ndm3=0;
    [ef,dis]=size(Red);
    Red_e=zeros(ef,dis);
    Red_m=zeros(ef,dis);
    Red_d=zeros(ef,dis);
    Red_em=zeros(ef,dis);
   
   [mic,dis]=size(Rmd);
    Rmd_e=zeros(mic,dis);
    Rmd_m=zeros(mic,dis);
    Rmd_d=zeros(mic,dis);
    Rmd_em=zeros(mic,dis);
	
	  [ef,mic]=size(Rem);
    Rem_e=zeros(ef,mic);
    Rem_m=zeros(ef,mic);
    Rem_d=zeros(ef,mic);
    Rem_dm=zeros(ef,mic);

        if(t<=l2)
        Red_e=a*E_1*Red+(1-a)*Y_2;   % ef -disease
        Red_m= Y_3 * Rmd;                 % ef-micrna *micRNA-disease
        nm2=1;
        ne2=1;
       end
    if(t<=r2)
        Red_d=a*Red*D_1+(1-a)*Y_2;   % ef-disease
		Red_em=Rem*Y_1;   % ef-mincrRNA* MicroRNA-disease
        nd2=1;
		nem2=1;
    end
    
    if(nm2+nd2+ne2+nem2>0)
    p4=ne2/(nm2+nd2+ne2+nem2);
    p5=nm2/(nm2+nd2+ne2+nem2);
    p6=nd2/(nm2+nd2+ne2+nem2);
	p7=nem2/(nm2+nd2+ne2+nem2);
     Red=p4*Red_e+p5*Red_m+p6*Red_d+p7*Red_em;
     else
        p4=0;
        p5=0;
        p6=0; 
		p7=0; 
        Red=Red;
    end
   
      
    if(t<=l1)
        Rmd_m=a*M_1*Rmd+(1-a)*Y_1;   % micrna -disease
        Rmd_e=Y_31* Red;                 %micrna -ef *EF-disease
        nm1=1;
        if(l2>0)
        ne1=1;
        else
            ne1=0;
        end
       
    end
    if(t<=r1)
        Rmd_d=a*Rmd*D_1+(1-a)*Y_1;   % disease-micrna
		Rmd_em=Rem'*Y_2; %microrna-ef * ef-disease
        nd1=1;
		nem1=1;
    end
    
    if(nm1+nd1+ne1+nem1>0)
     p1=nm1/(nm1+nd1+ne1+nem1);
    p2=nd1/(nm1+nd1+ne1+nem1);
    p3=ne1/(nm1+nd1+ne1+nem1);
	p8=nem1/(nm1+nd1+ne1+nem1);
    Rmd=p1*Rmd_m+p2*Rmd_d+p3*Rmd_e+p8*Rmd_em;
    else
        p1=0;
        p2=0;
        p3=0;  
		p8=0;  
         Rmd=Rmd;
    end
   
  
  
   if(t<=l3)
        Rem_m=a*Rem*M_1+(1-a)*Y_3;   % ef-micrna 
        Rem_ed=Red*Y_1';                 %EF-disease *disease_micrna
        nm3=1;
        ned3=1;
       
    end
    if(t<=r3)
        Rem_e=a*E_1*Rem+(1-a)*Y_3;   % ef-micrna
		Rem_dm=Y_2*Rmd'; % ef-disease*dis_microrna 
        ne3=1;
		nd3=1;
    end
    
    if(nm3+nd3+ne3+ned3>0)
     p9=nm3/(nm3+nd3+ne3+ned3);
    p10=nd3/(nm3+nd3+ne3+ned3);
    p11=ne3/(nm3+nd3+ne3+ned3);
	p12=ned3/(nm3+nd3+ne3+ned3);
    Rem=p9*Rem_d+p12*Rem_ed+p10*Rem_dm+p11*Rem_e;
    else
        p9=0;
        p10=0;
        p11=0;  
		p12=0;  
         Rem=Rem;
    end
   
end
    state.Rmd=Rmd;
    state.Red=Red;
	state.Rem=Rem;
end

