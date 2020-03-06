//
//  osci_e.C
//
//  Created by Tian XIE on 2/26/20.
//
//  This is to try getting three plots of neutrino oscillation. Inverted order.

#include "Tmath.h"
#include "TComplex.h"
#include "TMatrix.h"

void osci_e(){

    Double_t L,E;
    L=810.;
    E=1.8;

    //constants from particle data book
    Double32_t s12, s23, s13, c12, c23, c13, delta_21_s, delta_32_s;
    s12 = pow(0.307,0.5); s23 = pow(0.421,0.5); s13 = pow(2.12*pow(10.,-2.),0.5);
    c12 = cos(asin(s12)); c23 = cos(asin(s23)); c13 = cos(asin(s13));
    delta_21_s = 7.53*pow(10.,-5.); delta_32_s = 2.51*pow(10.,-3.);
    TMatrixD U(3,3);
    U(0,0)=c12*c13;
    U(0,1)=s12*c13;
    U(0,2)=s13;
    U(1,0)=-s12*c23-c12*s23*s13;
    U(1,1)=c12*c23-s12*s23*s13;
    U(1,2)=s23*c13;
    U(2,0)=s12*s23-c12*c23*s13;
    U(2,1)=-c12*s23-s12*c23*s13;
    U(2,2)=c23*c13;

    TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();
    TCanvas *c3 = new TCanvas();

    //p-E
    TF1 *p1 = new TF1("p1","1-4*([0]^2)*pow([1],2)*pow(sin(1.27*[3]*[5]/x),2)-4*pow([1],2)*pow([2],2)*pow(sin(1.27*([4]-[3])*[5]/x),2)-4*pow([0],2)*pow([2],2)*pow(sin(1.27*[4]*[5]/x),2)",0,5);
    p1->SetParameter(0,U(1,0));
    p1->SetParameter(1,U(1,1));
    p1->SetParameter(2,U(1,2));
    p1->SetParameter(3,delta_21_s);
    p1->SetParameter(4,delta_32_s);
    p1->SetParameter(5,L);
    p1->SetTitle("Survival Probability for #nu_{#mu} versus E (L=810km)");
    p1->GetHistogram()->GetXaxis()->SetTitle("E(GeV)");
    p1->GetHistogram()->GetYaxis()->SetTitle("Probability");
    c1->cd();
    p1->Draw();
    c1->Modified();
    c1->Update();


    //p-L
    TF1 * p2 = new TF1("p2","1-4*pow([0],2)*pow([1],2)*pow(sin(1.27*[3]*x/[5]),2)-4*pow([1],2)*pow([2],2)*pow(sin(1.27*([4]-[3])*x/[5]),2)-4*pow([0],2)*pow([2],2)*pow(sin(1.27*[4]*x/[5]),2)",0,1000);
    p2->SetParameter(0,U(1,0));
    p2->SetParameter(1,U(1,1));
    p2->SetParameter(2,U(1,2));
    p2->SetParameter(3,delta_21_s);
    p2->SetParameter(4,delta_32_s);
    p2->SetParameter(5,E);
    p2->SetTitle("Survival Probability for #nu_{#mu} versus L (E=1.8GeV)");
    p2->GetHistogram()->GetXaxis()->SetTitle("L(km)");
    p2->GetHistogram()->GetYaxis()->SetTitle("Probability");
    c2->cd();
    p2->Draw();
    c2->Modified();
    c2->Update();

    //p-L/E
    TF1 * p3 = new TF1("p3","1-4*pow([0],2)*pow([1],2)*pow(sin(1.27*[3]*x),2)-4*pow([1],2)*pow([2],2)*pow(sin(1.27*([4]-[3])*x),2)-4*pow([0],2)*pow([2],2)*pow(sin(1.27*[4]*x),2)",0,30000);
    p3->SetParameter(0,U(1,0));
    p3->SetParameter(1,U(1,1));
    p3->SetParameter(2,U(1,2));
    p3->SetParameter(3,delta_21_s);
    p3->SetParameter(4,delta_32_s);
    p3->SetTitle("Survival Probability for #nu_{#mu} versus L/E");
    p3->GetHistogram()->GetXaxis()->SetTitle("L/E(km/GeV)");
    p3->GetHistogram()->GetYaxis()->SetTitle("Probability");
    c3->cd();
    p3->Draw();
    c3->Modified();
    c3->Update();
    
}
