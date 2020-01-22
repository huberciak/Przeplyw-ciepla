#include<bits/stdc++.h>
#include <stdlib.h>
#include <algorithm>
using namespace std;
const int maxInt = INT_MAX - 1;
double fi (double x, double y)
{
   if (abs(x+1.0)<0.001) return -y;
   if (abs(x-1.0)<0.001) return y;
   if (abs(y-1.0)<0.001) return x;
   if (abs(y+1)<0.001) return -x;
}

void GaussElimination(vector<vector<double> >& resultMatrix, vector<double>& constant){

    double maxCol;
    int maxColNo;
    double proportion;

    for(int i = 0; i < resultMatrix.size(); i++){ //maksymalny element w kolumnie, kazdy wiersz
        maxColNo = i;
        maxCol = resultMatrix[i][i];

        for(int j = i+1; j<resultMatrix[i].size(); j++){
            if(resultMatrix[i][j] > maxCol){
                maxColNo = j;
                maxCol = resultMatrix[i][j];
            }
        }
        for(int j = i; j < resultMatrix.size(); j++){           //wiersz maksymalny idzie na i-te miejsce
            swap(resultMatrix[j][i], resultMatrix[j][maxColNo]);
        }
        swap(constant[i], constant[maxColNo]);

        if(resultMatrix[i][i] != 0){
            for(int j = i+1;j<resultMatrix.size();j++){          //wiersze wszystkie
                proportion = (resultMatrix[i][j])/(resultMatrix[i][i]);

                for(int k = i+1; k<resultMatrix.size(); k++){          //elementy wszystkie
                    resultMatrix[k][j] == resultMatrix[k][j] - resultMatrix[k][i]*proportion;
                }
                resultMatrix[i][j] = 0;
                constant[j] = constant[j] - constant[i]*proportion;
            }
        }                                       //MACIERZ Schodkowa
    }
}

vector< vector<double> > calculateMatrix (int size_) //glowna macierz do obliczen, Algorytm sekwencyjny generacji uk³adu równañ, wyliczenie B
{
//parametry startowe
    int pointsNumber = 3*size_*size_+4*size_+1;
    int shift =(2*size_+1) * (size_+1);//cout<<"shift = "<<shift<<endl;
    double a = 2.0/3, b = (-1.0/6), c = (-1.0/3); //jakies losowe calki
    double calculationMatrix [4][4] = { {a,b,c,b},{b,a,b,c},{c,b,a,b},{b,c,b,a}  };
//parametry startowe
//
//zerowanie macierzy
    vector< vector<double> > resultMatrix;
    resultMatrix.resize(pointsNumber);
    for(int i = 0; i < pointsNumber; i++)
    {
        resultMatrix[i].resize(pointsNumber);
        for (int j = 0; j < pointsNumber; j++) resultMatrix[i][j] = 0;

    }
//zerowanie macierzy
//
//Algorytm sekwencyjny generacji uk³adu równañ
    for (int i = 0; i <2*size_; i++)
    {
        for (int j = 0; j < size_; j++)
        { //cout<<"j ="<<j<<"i ="<<i<<"wykonianie"<<endl;
            int t[4] = {(j+1)*(2*size_+1)+i,(j+1)*(2*size_+1)+i+1,j*(2*size_+1)+i+1,j*(2*size_+1)+i}; //vector do odpowiedniego skakania po plytce

            for (int k = 0; k < 4; k++)
            {
                int Ek = pointsNumber - t[k] - 1;// cout<<"Ek ="<<Ek<<endl;
                resultMatrix[Ek][Ek] += calculationMatrix[k][k]; //cout<<"Ek = "<<Ek<<" Ek = "<<Ek<<" add = "<<calculationMatrix[k][k]<<" k = "<<k<<endl;
                for (int l = k+1; l<4; l++)
                {
                    int El = pointsNumber - t[l] - 1;
                resultMatrix[Ek][El] += calculationMatrix[k][l]; //cout<<"Ek = "<<Ek<<" El = "<<El<<" add = "<<calculationMatrix[k][l]<<endl;
                resultMatrix[El][Ek] += calculationMatrix[l][k]; //cout<<"El = "<<El<<" Ek = "<<Ek<<" add = "<<calculationMatrix[l][k]<<endl;
                }
                //cout<<"one round"<<endl;
            }
        }
    }
    //cout<<"CZESC DRUGA"<<endl;
    for (int i = 0; i < size_; i++)
    {
        for (int j = 0; j < size_; j++)
        {
            int t[4] = {(j*(size_+1)+i+shift),(j*(size_+1)+i+1+shift),((j-1)*(size_+1)+i+1+shift),((j-1)*(size_+1)+i+shift)};
            for (int k = 0; k < 4; k++)
            {
                int Ek = pointsNumber - t[k] - 1;
                resultMatrix[Ek][Ek] += calculationMatrix[k][k]; //cout<<"Ek = "<<Ek<<" Ek = "<<Ek<<" add = "<<calculationMatrix[k][k]<<" k = "<<k<<endl;
                for (int l = k+1; l < 4; l++)
                {
                    int El = pointsNumber - t[l] - 1;
                    resultMatrix[Ek][El] += calculationMatrix[k][l]; //cout<<"Ek = "<<Ek<<" El = "<<El<<" add = "<<calculationMatrix[k][l]<<endl;
                    resultMatrix[El][Ek] += calculationMatrix[l][k]; //cout<<"El = "<<El<<" Ek = "<<Ek<<" add = "<<calculationMatrix[l][k]<<endl;
                }
            }
        }
    }
//Algorytm sekwencyjny generacji uk³adu równañ
//
//zerowanie 4,5,7
    for(int i = size_*(2*size_+1); i < shift-size_; i++)
    {
        for(int j = 0; j < pointsNumber; j++)resultMatrix[j][i] = 0;
        resultMatrix[i][i] = 1; //cout<<"zerowane "<<i<<endl;
    }

    for (int i = (pointsNumber-1)/2 - 2*size_; i > 0; i -= (size_+1))
    {
        for (int j = 0; j < pointsNumber; j++)resultMatrix[j][i]=0;
        resultMatrix[i][i] = 1; //cout<<"zerowane "<<i<<endl;

    }
//zerowanie 4,5,7
//
//
    return resultMatrix;
}

vector<double> calculateL (int size_) //macierz wyrazow wolnych, wyliczamy L ;;; Algorytm sekwencyjny generacji uk³adu równañ czesc zwiazana z L ;;; fi zalezne od wiersza
{
    int p = 3*size_*size_+4*size_ + 1;
    double div = 1.0/size_;
    int n = 2*size_;
    vector <double> L;
    L.resize(p);

    for(int i = 0; i < p; i++) L[i] = 0;

    L[0] = 0.5* fi( -1 , (1-(0.5)*div) )*div + 0.5*fi( (-1+(0.5)*div) , 1 )*div;                       //warunki brzegowe
    L[n] = 0.5* fi( 1 , (1-0.5*div) )*div + 0.5*fi( (1-0.5*div) , 1 )*div;
    L[p - 1] = 0.5* fi( 1 , (-1+0.5*div) )*div + 0.5*fi( (1-0.5*div) , -1 )*div;

	for(int i=1; i < n; i++) L[i] = 0.5* fi( (-1-0.5*div+i*div), 1)*div + 0.5*fi( (-1+0.5*div+i*div), 1)*div;

    n++;

	for(int i = 1; i < size_; i++)
    {
        L[i*n] = 0.5* fi( -1 , (1+0.5*div-i*div) )*div + 0.5* fi( -1 , (1-0.5*div-i*div) )*div;
        L[i*n+2*size_] = 0.5* fi( 1 , (1+0.5*div-i*div) )*div + 0.5* fi( 1 , (1-0.5*div-i*div) )*div;
    }

    n *= (size_+1) - 1;

    for(int i = 0; i < size_; i++) L[n+i*(size_+1)] = 0.5*fi( 1 , (0.5*div-i*div) )*div + 0.5*fi( 1 , (-0.5*div-i*div) )*div;

    for (int i = 0; i < size_-1; i++) L[p-size_+i] = 0.5*fi( (0.5*div+i*div) , -1 )*div + 0.5*fi( (1.5*div+i*div) , -1 )*div;

    return L;
}

vector<double> calculateVectorW(vector<vector<double> > B, vector<double> L){

    vector<double> resultVector;
    resultVector.resize(B.size());

    for(int i = B.size()-1; i >= 0; i--){                    //idac od ostatniego wiersza odejmujemy wszystkie poprzednie wartosci bo jest zgaussowane
        resultVector[i] = L[i];
        for(int j = B.size()-1; j > i; j--){
            resultVector[i] -= B[j][i]*resultVector[j];
        }

        if(B[i][i] != 0) resultVector[i] = resultVector[i]/B[i][i];
        else
            resultVector[i] = maxInt;
    }
    return resultVector;
}

int main()
{
    int N;        //ilosc podzialow
    cin >> N;
    int points = 3*N*N +4*N + 1;

    vector < vector<double> > B = calculateMatrix(N);
    vector <double> L = calculateL(N);
    GaussElimination(B, L);
    vector <double> result;
    result = calculateVectorW(B, L);

    cout<<"[ ";
    for (int i=0; i<result.size(); i++)
    {
        if(i!=result.size()-1)
        cout<<result[i]<<", ";
        else
            cout<<result[i]<<" ";         //wektor w
    }
    cout<<"]"<<endl;;





    vector< vector<double> > testresult = calculateMatrix(N);
    for(int i=0;i<points;i++){
        cout<<endl;
        for(int j=0;j<points;j++){
                cout<<testresult[i][j]<<" ";

        }

    }

}





