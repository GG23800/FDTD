#ifndef HEAT_FDTD_1D_HPP
#define HEAT_FDTD_1D_HPP

#include<string>
#include<chrono>
#include<iostream>
#include<fstream>
#include<thread>
#include<mutex>

#include<unistd.h>

#include"Vector.hpp"
#include"material.hpp"

/*
 * Contributor: Jérôme Dubois
 * Date: 06/2022
 * Version: 0.1
 *
 *
 */

enum BoundaryCondition
{
    None = 0,
    Dirichlet,
    Neumann,
    Periodic,
    NbBC
};

enum SaveFileType
{
    txt = 0,
    json,
    binary,
    Nb
};
// TODO add FDTDError ???
/*
enum FDTDError
{
    NoError = 0;
    BoundaryCondition;
    MaterialIDExceeding;
    DeltatToHigh

};
*/

// TODO gerer la valeur et l'ajout de materiel
// TODO reset loopcount quand changement

class FDTD_H1D
{
    public:
        MaterialList<HeatMaterial> HMList;

        FDTD_H1D();
        ~FDTD_H1D();

        bool IsValidState() {return Valid;}

        // sources manipulation
        void ResetSources(float value = 0.f);
        void SetSources(unsigned int node, float value);
        void SetSources(unsigned int node0, unsigned int nodef, float value);
        void SetSources(Vector<float> nsources);
        Vector<float> GetSources() {return sources;}

        // Material manipulation
        void ResetMaterials(unsigned short mat = 0);
        void SetMaterials(unsigned int node, unsigned short mat);
        void SetMaterials(unsigned int node0, unsigned int nodef, unsigned short mat);
        void SetMaterials(Vector<unsigned short> nmaterials);
        Vector<unsigned short> GetMaterials() {return materials;}

        // Initial condition
        void ResetInitialCondition(float value = 0.f);
        void SetRandomInitialCondition(float Tmin, float Tmax);
        void SetInitialCondition(unsigned int node, float value);
        void SetInitialCondition(unsigned int node0, unsigned int nodef, float value);
        void SetInitialCondition(Vector<float> nic);
        Vector<float> GetInitialCondition() {return evenVector;}

        void SetPeriodicBoundaryCondition();
        void SetBoundaryCondition(unsigned int node, BoundaryCondition nbc, float value = 0.f);

        // Convergence test function
        void SetAutoConvergenceTestOn() {AutoConvergenceTest = true;}
        void SetAutoConvergenceTestOff() {AutoConvergenceTest = false;}
        bool GetAutoConvergenceTestState() {return AutoConvergenceTest;}
        void SetConvergenceTestStep(unsigned int ncts) {ConvergenceTestStep = ncts;}
        unsigned int GetConvergenceTestStep() {return ConvergenceTestStep;}
        void SetConvergenceTestCondition(float nctc);
        float GetConvergenceTestCondition() {return ConvergenceTestCondition;}
        bool GetConvergenceResult() {return Converged;}

        // Periodic save functions
        void SetPeriodicSaveFileType(SaveFileType nSaveFileType);
        void SetPeriodicSaveStep(unsigned int nSaveStep) {SaveStep =nSaveStep;}
        void SetPeriodicSaveState(bool nPeriodicSaveState) {PeriodicSave = nPeriodicSaveState;}
        bool GetPeriodicSaveState() {return PeriodicSave;}
        void PeridodicSaveStep(unsigned int nSaveStep);
        void SetPeriodicSaveFileName(std::string nName) {PeriodicSaveFileName = nName; PeriodicSaveFileRenamed = true;}
        
        void SetDeltax(float dx);
        void SetDeltat(float dt);
        void SetNumberOfNode(unsigned int nnon);
        void SetNumberOfCalculationLoop(unsigned int nnocl);
        void ResetValidity() {Valid = true;} // This function is dedicated to reset valid state to true (TODO add automatic reset on modification?). This function is intented to be called before changing simulation settings (in case one simulation has runned correctly and user want to iterate on a new simulation for convergence test for example, in this case actually user HAVE to use this function). Actually it can't be called to force simulation validity, a sanity check is called automatically before running a simulation

        float GetDeltax() {return Deltax;}
        float GetDeltat() {return Deltat;}
        unsigned int GetNumberOfNode() {return NumberOfNode;}
        unsigned int GetNumberOfCalculationLoop() {return MaxLoop;}

        void CheckForError(); // TODO return last error?
        void Run();
        Vector<float> GetLastSimulationOutput();
        void GetLastSimulationOutput(Vector<float> Output);
        Vector<float> GetNodeWeight() {return NodeWeight;}
        inline void Stop() {Running = false;}
        inline bool IsRunning() {return Running;}
        inline bool IsPreparing() {return Preparing;}

        // TODO check list
        // implement dynamic sources (and later/on the same time any other dynamic parameter)
        // implement a get last simulation result, need mutex
        // implement a periodic

    private:
        void CheckDeltat();
        void CheckBoundaryCondition();
        void CheckMaterialValidity();
        void CheckMediumThickness();

        void PrepareBoundaryCalculation();
        void PrepareNodeWeightForDeltat();
        void PrepareNodeWeight();
        void AdaptInterfaceWeight();
        // TODO ? void AdaptSources();
        void PrepareKernel();

        void CalculateBoundary(); // TODO check function...
        void prun();

        // autosave
        void OpenPeriodicFile();
        void ClosePeriodicFile();
        void WriteActualStep();
        std::string GetAutosaveFileName();

        float Deltax;
        float Deltat;
        unsigned int NumberOfNode;

        BoundaryCondition LeftBC;
        BoundaryCondition RightBC;
        float LeftBCValue;
        float RightBCValue;
        Vector<unsigned int> LeftID;
        Vector<unsigned int> RightID;
        Vector<float> LeftWeight; // TODO remove? pas utile?
        Vector<float> RightWeight; // TODO idem

        unsigned long int MaxLoop;
        unsigned long int LoopCount;
        unsigned long int LastAskedOutput;
        bool Dynamic;
        bool Valid;
        bool Preparing;
        bool Running;
        std::recursive_mutex OddVectorProtection;
        std::recursive_mutex EvenVectorProtection;

        // Autotest
        bool AutoConvergenceTest;
        bool Converged;
        unsigned int ConvergenceTestStep;
        float ConvergenceTestCondition;

        // Save options
        // save state in json
        // save type for periodic save
        bool PeriodicSave;
        bool PeriodicSaveFileRenamed;
        SaveFileType lSaveFileType;
        unsigned int SaveStep;
        std::string PeriodicSaveFileName;
        std::string PeriodicSaveExtension;
        std::ofstream PeriodicSaveFile;

        Vector<float> sources;
        Vector<unsigned short> materials;
        Vector<float> evenVector;
        Vector<float> oddVector;
        Vector<float> NodeWeight;
        Vector<float> *ActiveVector;
        Vector<float> *RefVector;
};

#endif
