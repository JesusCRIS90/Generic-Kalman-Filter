#ifndef KALMANFILTER_HPP
#define KALMANFILTER_HPP

#include <armadillo/include/armadillo>
#include <iostream>

typedef arma::mat BasicMatrix;
typedef arma::vec BasicVector;

//enum VerifyType { VerifyVector, Verify_A, Verify_B, Verify_C, Verify_D, Verify_Q, Verify_R, Verify_P};



class KalmanVerify{

public:
    static void Verify(BasicMatrix&, arma::uword, arma::uword);
    static void Verify(BasicVector&, arma::uword);
private:
    KalmanVerify();
};



class KalmanParameters
{
public:
    KalmanParameters(BasicVector, BasicVector, BasicMatrix, BasicMatrix, BasicVector = {}, BasicMatrix = {});
    KalmanParameters(BasicMatrix, BasicMatrix, BasicMatrix, BasicMatrix, BasicMatrix,
                     BasicVector, BasicVector, BasicMatrix = {}, BasicMatrix = {}, BasicVector = {});
    KalmanParameters(KalmanParameters& ptr);
    KalmanParameters(KalmanParameters* ptr);
    KalmanParameters(KalmanParameters&& ptr);

    KalmanParameters& set_A(BasicMatrix = {});
    KalmanParameters& set_B(BasicMatrix = {});
    KalmanParameters& set_C(BasicMatrix = {});
    KalmanParameters& set_D(BasicMatrix = {});
    KalmanParameters& set_Q(BasicMatrix = {});
    KalmanParameters& set_R(BasicMatrix = {});

    KalmanParameters& set_Pk(BasicMatrix = {});

    void copyParameters(KalmanParameters* data2copy);
    void copyParameters(KalmanParameters& data2copy);


    friend class KalmanFilter;

    friend std::ostream& operator << (std::ostream&, const KalmanParameters&);
    friend std::ostream& operator << (std::ostream&, const KalmanFilter&);
//    const KalmanParameters& operator= (KalmanParameters& Obj2Copy);

private:
    // constant over all executation
    BasicMatrix A = {};
    BasicMatrix B = {};
    BasicMatrix C = {};
    BasicMatrix D = {};
    BasicMatrix Q = {};
    BasicMatrix R = {};

    BasicMatrix Pk_k_corregida = {};             // Need a Set
    BasicVector Xk_k_corregido = {};             // It is the OUTPUT

    // variable over all executation
    BasicMatrix Pk_km1_estimada  = {};           // Internal Matrix
    BasicMatrix KalmanGain = {};                 // Internal Matrix
    BasicVector Xk_km1_estimado = {};            // Internal Vector
    BasicVector Uk = {};                         // It is neccesary
    BasicVector Yk_estimada = {};                // It is the INPUT


    arma::uword Nx = 0;             // Size of Vectors and Matrix. Most important parameter
    arma::uword Nz = 0;
    arma::uword Nu = 0;
};


class KalmanFilter
{
public:
    //KalmanFilter();
    KalmanFilter(BasicVector, BasicVector, BasicMatrix, BasicMatrix, BasicVector = {}, BasicMatrix = {});
    KalmanFilter(BasicMatrix, BasicMatrix, BasicMatrix, BasicMatrix, BasicMatrix,
                 BasicVector, BasicVector, BasicMatrix = {}, BasicMatrix = {}, BasicVector = {});

    KalmanFilter(KalmanParameters);
    ~KalmanFilter();

    void Init(BasicVector, BasicVector, BasicMatrix, BasicMatrix, BasicVector = {}, BasicMatrix = {});
    void Init(BasicMatrix, BasicMatrix, BasicMatrix, BasicMatrix, BasicMatrix,
              BasicVector, BasicVector, BasicMatrix = {}, BasicMatrix = {}, BasicVector = {});

    void update(BasicVector input, BasicVector = {} );      // Most Important function. Calling each iteration.

    BasicVector GetState(){ return params->Xk_k_corregido;}
    BasicVector GetMeasure(){ return params->Yk_estimada;}

    void SetConfiguration(KalmanParameters);
    KalmanParameters Configuration();

private:

    void CorrectionPhase();
    void PredictionPhase();

    //friend class KalmanParameters;
    friend std::ostream& operator<< (std::ostream&, const KalmanFilter&);

    KalmanParameters* params;
    BasicVector current_input = {};

};

#endif // KALMANFILTER_HPP
