#include"heat_fdtd_1d.hpp"

FDTD_H1D::FDTD_H1D() : LeftID(0), RightID(0), LeftWeight(0), RightWeight(0), sources(3), materials(3), evenVector(3), oddVector(3), NodeWeight(3)
{
    Deltax = 1.f;
    Deltat = 1.f;
    NumberOfNode = 3;

    LeftBC = Dirichlet;
    RightBC = Dirichlet;
    LeftBCValue = 0.f;
    RightBCValue = 0.f;

    MaxLoop = 10;
    LoopCount = 0;
    Dynamic = false;
    Valid = true;
    Running = false;

    // Convergence test
    AutoConvergenceTest = false;
    Converged = true;
    ConvergenceTestStep = 1000;
    ConvergenceTestCondition = 0.001f;

    sources.reset();
    materials.reset();
    evenVector.reset();
    oddVector.reset();
    NodeWeight.reset();
    ActiveVector = &evenVector;
    RefVector = &oddVector;
}

FDTD_H1D::~FDTD_H1D()
{
    if (Running)
    {
        Running=false;
        // TODO wait or add a bool to be sure that if a calculation loop
        // we wait it is finished or a sleeping time
    }
}

void FDTD_H1D::ResetSources(float value)
{
    sources.reset(value);
}

void FDTD_H1D::SetSources(unsigned int node, float value)
{
    if (node >= NumberOfNode)
    {
        std::cout << "Warning, trying to change sources at node " << node << " while there are only " << NumberOfNode << " nodes... Nothing is done." << std::endl;
    }
    else
    {
        sources.data[node] = value;
    }
}

void FDTD_H1D::SetSources(unsigned int node0, unsigned int nodef, float value)
{
    if (nodef < node0)
    {
        unsigned int tmp = nodef;
        nodef = node0;
        node0 = tmp;
    }
    if (node0 >= NumberOfNode) 
    {
        std::cout << "Warning, trying to change sources from node " << node0 << " to node " << nodef << " while there are only " << NumberOfNode << " nodes... Nothing is done." << std::endl;
    }
    else if (nodef >= NumberOfNode)
    {
         std::cout << "Warning, trying to change sources from node " << node0 << " to node " << nodef << " while there are only " << NumberOfNode << " nodes... Sources value is changed only for valid nodes." << std::endl;
         for (unsigned int k=node0 ; k<NumberOfNode ; k++)
         {
            sources.data[k] = value;
         }
    }
    else
    {
         for (unsigned int k=node0 ; k<nodef ; k++)
         {
            sources.data[k] = value;
         }       
    }
}

void FDTD_H1D::SetSources(Vector<float> nsources)
{
    if (nsources.size() != NumberOfNode)
    {
        std::cout << "Warning, trying to change sources with a vector of length " << nsources.size() << " whereas the simulation has " << NumberOfNode << ". Nothing is done." << std::endl;
    }
    else
    {
        sources = nsources;
    }
}

void FDTD_H1D::ResetMaterials(unsigned short mat)
{
    materials.reset(mat);
}

void FDTD_H1D::SetMaterials(unsigned int node, unsigned short mat)
{
    if (node >= NumberOfNode)
    {
        std::cout << "Warning, trying to change materials at node " << node << " while there are only " << NumberOfNode << " nodes... Nothing is done." << std::endl;
    }
    else
    {
        materials.data[node] = mat;
    }
}

void FDTD_H1D::SetMaterials(unsigned int node0, unsigned int nodef, unsigned short mat)
{
    if (nodef < node0)
    {
        unsigned int tmp = nodef;
        nodef = node0;
        node0 = tmp;
    }
    if (node0 >= NumberOfNode) 
    {
        std::cout << "Warning, trying to change materials from node " << node0 << " to node " << nodef << " while there are only " << NumberOfNode << " nodes... Nothing is done." << std::endl;
    }
    else if (nodef >= NumberOfNode)
    {
         std::cout << "Warning, trying to change materials from node " << node0 << " to node " << nodef << " while there are only " << NumberOfNode << " nodes... Sources value is changed only for valid nodes." << std::endl;
         for (unsigned int k=node0 ; k<NumberOfNode ; k++)
         {
            materials.data[k] = mat;
         }
    }
    else
    {
         for (unsigned int k=node0 ; k<nodef ; k++)
         {
            materials.data[k] = mat;
         }       
    }
}

void FDTD_H1D::SetMaterials(Vector<unsigned short> nmaterials)
{
    if (nmaterials.size() != NumberOfNode)
    {
        std::cout << "Warning, trying to change materials with a vector of length " << nmaterials.size() << " whereas the simulation has " << NumberOfNode << ". Nothing is done." << std::endl;
    }
    else
    {
        materials = nmaterials;
    }
}

void FDTD_H1D::ResetInitialCondition(float value)
{
    evenVector.reset(value);
}

void FDTD_H1D::SetRandomInitialCondition(float Tmin, float Tmax)
{
    evenVector.SetRandom(Tmin, Tmax);
}

void FDTD_H1D::SetInitialCondition(unsigned int node, float value)
{
    if (node >= NumberOfNode)
    {
        std::cout << "Warning, trying to change initial condition at node " << node << " while there are only " << NumberOfNode << " nodes... Nothing is done." << std::endl;
    }
    else
    {
        evenVector.data[node] = value;
    }
}

void FDTD_H1D::SetInitialCondition(unsigned int node0, unsigned int nodef, float value)
{
    if (nodef < node0)
    {
        unsigned int tmp = nodef;
        nodef = node0;
        node0 = tmp;
    }
    if (node0 >= NumberOfNode) 
    {
        std::cout << "Warning, trying to change initial condition from node " << node0 << " to node " << nodef << " while there are only " << NumberOfNode << " nodes... Nothing is done." << std::endl;
    }
    else if (nodef >= NumberOfNode)
    {
         std::cout << "Warning, trying to change initial condition from node " << node0 << " to node " << nodef << " while there are only " << NumberOfNode << " nodes... Sources value is changed only for valid nodes." << std::endl;
         for (unsigned int k=node0 ; k<NumberOfNode ; k++)
         {
            evenVector.data[k] = value;
         }
    }
    else
    {
         for (unsigned int k=node0 ; k<nodef ; k++)
         {
            evenVector.data[k] = value;
         }       
    }
}

void FDTD_H1D::SetInitialCondition(Vector<float> nic)
{
    if (nic.size() != NumberOfNode)
    {
        std::cout << "Warning, trying to change sources with a vector of length " << nic.size() << " whereas the simulation has " << NumberOfNode << ". Nothing is done." << std::endl;
    }
    else
    {
        evenVector = nic;
    }
}

void FDTD_H1D::SetPeriodicBoundaryCondition()
{
    LeftBC = Periodic;
    RightBC = Periodic;
}

void FDTD_H1D::SetBoundaryCondition(unsigned int node, BoundaryCondition nbc, float value)
{
    if (node >= NumberOfNode)
    {
        std::cout << "Warning, trying to edit boundary condition for an out of range node, right boundary condition will be edited." << std::endl;
        node = NumberOfNode-1;
    }
    if (node == 0)
    {
        LeftBC = nbc;
        LeftBCValue = value;
    }
    else if (node == NumberOfNode-1)
    {
        RightBC = nbc;
        RightBCValue = value;
    }
    else
    {
        std::cout << "Warning, internal boundary condition not implemented yet..." << std::endl;
    }
}

void FDTD_H1D::SetConvergenceTestCondition(float nctc)
{
    if (nctc <= 0.f)
    {
        std::cout << "Warning, the convergence condition value can't be less or equal to 0. It is set to it's default value: 1e-3. This may lead to error later." << std::endl;
        ConvergenceTestCondition = 1e-3f;
    }
    else if (nctc < 1e-10f)
    {
        std::cout << "Warning, the convergence condition value of " << nctc << " is very small. Convergence can be very long or never happen, it is at your own risk." << std::endl;
    }
    else
    {
        ConvergenceTestCondition = nctc;
    }
}

void FDTD_H1D::CheckForError()
{
    CheckBoundaryCondition();
    CheckMaterialValidity();
    CheckMediumThickness();
    CheckDeltat();
}

void FDTD_H1D::SetDeltax(float dx)
{
    if (dx <= 0.f)
    {
        std::cout << "Warning, asked deltax " << dx << " can not be set, nothing is done..." << std::endl;
    }
    else
    {
        Deltax = dx;
    }
}

void FDTD_H1D::SetDeltat(float dt)
{
    if (dt <= 0.f)
    {
        std::cout << "Warning, asked deltat " << dt << " can not be set, nothing is done..." << std::endl;
    }
    else
    {
        Deltat = dt;
    }
}

void FDTD_H1D::SetNumberOfNode(unsigned int nnon)
{
    if (nnon < 2)
    {
        std::cout << "Warning, asked number of node " << nnon << " can not be set, nothing is done..." << std::endl;
    }
    else
    {
        NumberOfNode = nnon;
        // TODO change size of all dynamic vector
        sources.resize(nnon);
        sources.reset();
        materials.resize(nnon);
        materials.reset();
        evenVector.resize(nnon);
        evenVector.reset();
        oddVector.resize(nnon);
        oddVector.reset();
        NodeWeight.resize(nnon);
        NodeWeight.reset();
    }
}

void FDTD_H1D::SetNumberOfCalculationLoop(unsigned int nnocl)
{
    if (nnocl <1)
    {
        std::cout << "Warning, the number of calculation loop is set to 0. Automaticaly change to 1 to have at least one calculation." << std::endl;
        nnocl = 1;
    }
    MaxLoop = nnocl;
}


void FDTD_H1D::CheckDeltat()
{
    PrepareNodeWeightForDeltat();
    std::cout << "NodeWeight: ";
    NodeWeight.print();
    std::cout << "NodeWeight max: " << NodeWeight.max() << std::endl;
    float dtm = 0.5f/NodeWeight.max();
    std::cout << "Deltat max: " << dtm << std::endl;
    if (Deltat >= dtm) 
    {
        std::cout << "Warning, deltat is too high (" << Deltat << ") and the simulation will diverge. Deltat maximum value (not included) is " << dtm <<", we set it's value to " << 0.99f*dtm << ". This may lead to an error later." << std::endl;
        Deltat = 0.99f*dtm;
    }
    Valid = true;
    // TODO how to manage Valid?? Because if false, it can't become true in function becaus it can be false before entering the function
    // Moreover, actually this function will always lead to a true value of Valid parameter
    // Do we allways start with a true Valid boolean and put it to false only if valid and display the warning? This option will need a reset function at least for Valid parameter
    // Other possibility will be to make an enum on validity state with error, unknow and true/valid/no_error possibility.
    // In first implementation, CheckDeltat will allways be called first and whatever the state of Valid it will set it to true. The other Check funtion can only set Valid to false, not to true
}
// TODO implement a verbose option
// so here we will have a verbose condition
// can be a compile option
void FDTD_H1D::CheckBoundaryCondition()
{
    if (LeftBC == None) 
    {
        std::cout << "Error, no boundary condition affected to left node (ID=0)." << std::endl;
        Valid = false;
    }
    else if (RightBC == None) 
    {
        std::cout << "Error, no boundary condition affected to right node (ID=" << NumberOfNode-1 <<")." << std::endl;
        Valid = false;
    }
    else
    {
        if (LeftBC == Periodic)
        {
            if (RightBC == LeftBC) {}//Valid = true;}
            else 
            {
                std::cout << "Error, right boundary condition (node " << NumberOfNode-1 << ") is not periodic whereas left one is. Please verifiy your simulation definition." << std::endl;
                Valid = false;
            }
        }
        if (RightBC == Periodic)
        {
            if (LeftBC == RightBC) {}//Valid = true;}
            else 
            {
                std::cout << "Error, left boundary condition (node 0) is not periodic whereas right one is. Please verifiy your simulation definition." << std::endl;
                Valid = false;
            }
        }
    }
    // TODO error definition?? Like ErrorBoundaryDefinition?
}

void FDTD_H1D::PrepareBoundaryCalculation()
{
    // TODO finish with material values
    if (LeftBC == Dirichlet)
    {
        LeftID.resize(0);
        LeftWeight.resize(0);
    }
    else if (LeftBC == Neumann)
    {
        LeftID.resize(2);
        LeftWeight.resize(2);
        LeftID.data[0] = 0;
        LeftID.data[1] = 1;
    }
    else if (LeftBC == Periodic)
    {
        LeftID.resize(3);
        LeftWeight.resize(3);
        LeftBCValue = 0.f;
        RightBCValue = 0.f;
        LeftID.data[0] = NumberOfNode-2;
        LeftID.data[1] = 0;
        LeftID.data[2] = 1;
    }
    else
    {
        std::cout << "Warning unknow left Boundary Condition in PrepareBoundaryCalculation..." << std::endl;
    }

    if (RightBC == Dirichlet)
    {
        RightID.resize(0);
        RightWeight.resize(0);
    }
    else if (RightBC == Neumann)
    {
        RightID.resize(2);
        RightWeight.resize(2);
        RightID.data[0] = NumberOfNode-1;
    }
    else if (RightBC == Periodic)
    {
        RightID.resize(0);
        RightWeight.resize(0);
        //RightID.data[0] = NumberOfNode-1;
        //RightWeight.data[0] = 1.f;
    }
    else
    {
        std::cout << "Warning unknow right Boundary Condition in PrepareBoundaryCalculation..." << std::endl;
    }
}

void FDTD_H1D::CalculateBoundary()
{
    unsigned int k=0;
    unsigned int pos=0;
    if ( (LoopCount%2) == 0 ) // TODO only once at master calculation function
    {
        ActiveVector = &oddVector;
        RefVector = &evenVector;
    }
    else
    {
        ActiveVector = &evenVector;
        RefVector = &oddVector;
    }

    ActiveVector->data[pos] = 0.f;
    for (k=0 ; k<LeftID.size() ; k++)
    {
        ActiveVector->data[pos] += LeftWeight.data[k]*RefVector->data[LeftID.data[k]];
    }
    ActiveVector->data[pos] += LeftBCValue;

    pos = NumberOfNode-1;
    ActiveVector->data[pos] = 0.f;
    for (k=0 ; k<RightID.size() ; k++)
    {
        ActiveVector->data[pos] += RightWeight.data[k]*RefVector->data[RightID.data[k]];
    }
    ActiveVector->data[pos] += RightBCValue;
}

void FDTD_H1D::CheckMaterialValidity()
{
    short MaxMaterialID = materials.max();
    if ( (unsigned int)(MaxMaterialID) > HMList.size() )
    {
        std::cout << "Error at least one material assigned in the simulation is not defined in the material list (material ID > material list size)... Invalid state." << std::endl;
        Valid = false;
        // TODO LastError?
    }
    else
    {
        //Valid = true;
        // can't change valid to true. Because what if it was false before?
    }
}

void FDTD_H1D::CheckMediumThickness()
{
    // TODO error definition? Actually will not lead to invalid state, just launch a warning
    unsigned short LastInterface=materials.data[0];
    unsigned short NewInterface=materials.data[0];
    for (unsigned int k=1 ; k<NumberOfNode-1 ; k++)
    {
        if (materials.data[k] != materials.data[k-1])
        {
            if (k==1)
            {
                std::cout << "Warning, there is an interface at node 1, juste next to the left boundary condition, take at least 2 pixels. Simulation will run but is not well defined. An error may occur later." << std::endl;
            }
            NewInterface = k;
            if ( (NewInterface-LastInterface)==1 )
            {
                // TODO errorthickness?
                // TODO error if thickness <= 2 pixels or minimal thickness can be edited.
                std::cout << "Warning, there is a material plate/slice of only one pixel at node " << k << ". Simulation will run but it is not well defined. An error may occur later." << std::endl;
            }
            LastInterface = NewInterface;
        }
    }
    if (materials.data[NumberOfNode-1] != materials.data[NumberOfNode-2])
    {
        std::cout << "Warning, there is an interface at node " << NumberOfNode-1 << ", juste next to the right boundary condition, take at least 2 pixels. Simulation will run but is not well defined. An error may occur later." << std::endl;
    }
}

void FDTD_H1D::PrepareNodeWeightForDeltat()
{
    if (Valid)
    {
        float lThermDiff = 0.f;
        unsigned int lmID = 0;
        lmID = (unsigned int)(materials.data[0]);
        lThermDiff = (HMList.get_material(lmID)).get_thermal_diffusivity();
        NodeWeight.data[0] = lThermDiff/Deltax/Deltax;
        for (unsigned int k=1 ; k<NumberOfNode ; k++)
        {
            if ( (unsigned int)(materials.data[k] != lmID) )
            {
                lmID = (unsigned int)(materials.data[k]);
                lThermDiff = (HMList.get_material(lmID)).get_thermal_diffusivity();
                NodeWeight.data[k] = lThermDiff/Deltax/Deltax;
            }
            else
            {
                NodeWeight.data[k] = NodeWeight.data[k-1];
            }
        }
    }
}

void FDTD_H1D::PrepareNodeWeight()
{
    if (Valid)
    {
        PrepareNodeWeightForDeltat();
        // TODO add a check to not recalculate that?
        for (unsigned int k=0 ; k<NumberOfNode ; k++)
        {
            NodeWeight.data[k] *= Deltat;
        }
    }
    else
    {
        NodeWeight.reset(0.f);
        std::cout << "Impossible to prepare nodes weight, invalide state..." << std::endl;
    }
}

void FDTD_H1D::AdaptInterfaceWeight()
{
    // for clarity it is not made at the same time than CheckMediumThickness
    // wich is dedicated to check for error
    // This function must be called after PrepareNodeWeight()
    Vector<float> lNW(NodeWeight);
    for (unsigned int k=1 ; k<NumberOfNode-1 ; k++)
    {
        if (materials.data[k] != materials.data[k-1])
        {
            NodeWeight.data[k] = (lNW.data[k-1] + lNW.data[k+1])/2.f;
        }
    }
}

void FDTD_H1D::PrepareKernel()
{
    Valid = true;
    CheckForError();
    PrepareBoundaryCalculation();
    PrepareNodeWeight();
    AdaptInterfaceWeight();
    LoopCount = 0;
}

void FDTD_H1D::Run()
{
    std::thread::id this_id = std::this_thread::get_id();
    std::cout << "Main thread id: " << this_id << std::endl;
    CheckForError();
    PrepareKernel();
    std::thread run_thread(&FDTD_H1D::prun, this);
    run_thread.detach();
}

// TODO inclue an auto convergence test with step defined by user
void FDTD_H1D::prun()
{
    // TODO make it as a private function
    // the public function will send this function in a different thread
    // TODO include Mutex protection
    std::thread::id this_id = std::this_thread::get_id();
    std::cout << "prun id: " << this_id << std::endl;
    Converged = false;
    unsigned int k = 0;
    unsigned int NextConvergenceTest = ConvergenceTestStep;
    Vector<float> ConvergenceDiff(0);
    if (Valid)
    {
        Running = true;
        while ( Running && (LoopCount < MaxLoop) )
        {
            if (LoopCount%100 == 0)
            {
                std::cout << "Running loop: " << LoopCount << std::endl;
            }
            if ( (LoopCount%2) == 0 )
            {
                OddVectorProtection.lock();
                ActiveVector = &oddVector;
                RefVector = &evenVector;
            }
            else
            {
                EvenVectorProtection.lock();
                ActiveVector = &evenVector;
                RefVector = &oddVector;
            }

            CalculateBoundary();
            for (k=1 ; k<NumberOfNode-1 ; k++)
            {
                ActiveVector->data[k] = RefVector->data[k] + NodeWeight.data[k-1]*RefVector->data[k-1] + NodeWeight.data[k+1]*RefVector->data[k+1] - 2.f*NodeWeight.data[k]*RefVector->data[k] + sources.data[k];
            }
            LoopCount++;
            // LoopCount has been incremented so inversion on the test compare to first one
            if ( (LoopCount%2) == 0 ) {EvenVectorProtection.unlock();}
            else {OddVectorProtection.unlock();}
            
            if (AutoConvergenceTest)
            {
                // Actually there is convergence test in addition to the maximum number of loop. Maximum number of loop can act as a security to not finish in an infinite loop for convergence
                if (LoopCount == NextConvergenceTest)
                {
                    ConvergenceDiff = (*ActiveVector) - (*RefVector);
                    std::cout << "Convergence test at loop " << LoopCount << ", convergence diff max: " << ConvergenceDiff.AbsMax() << std::endl;
                    if (ConvergenceDiff.AbsMax() < ConvergenceTestCondition)
                    {
                        Running = false;
                        Converged = true;
                    }
                    NextConvergenceTest += ConvergenceTestStep;
                }
            }
        }
        Running = false;
    }
    else
    {
        std::cout << "Error, trying to run a simulation whereas it is not valid. Abortion." << std::endl;
    }
}

Vector<float> FDTD_H1D::GetLastSimulationOutput()
{
    // not very efficient leads to creation and copy of a vector
    // prefer second version of the functin for time improvement
    Vector<float> output(0);
    if ( (LoopCount%2) == 0 ) {OddVectorProtection.lock();}
    else {EvenVectorProtection.lock();}
    output = *(ActiveVector);
    if ( (LoopCount%2) == 0 ) {OddVectorProtection.unlock();}
    else {EvenVectorProtection.unlock();}
    return output;
}

void FDTD_H1D::GetLastSimulationOutput(Vector<float> Output)
{
    if ( (LoopCount%2) == 0 ) {OddVectorProtection.lock();}
    else {EvenVectorProtection.lock();}
    Output = *(ActiveVector);
    if ( (LoopCount%2) == 0 ) {OddVectorProtection.unlock();}
    else {EvenVectorProtection.unlock();}
}
