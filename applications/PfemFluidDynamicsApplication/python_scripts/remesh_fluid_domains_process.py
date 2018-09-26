from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
KratosMultiphysics.CheckForPreviousImport()


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshFluidDomainsProcess(Model, settings["Parameters"])


class RemeshFluidDomainsProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        self.main_model_part = Model[custom_settings["model_part_name"].GetString()]

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"            : 0,
            "model_part_name"       : "Fluid Domain",
            "meshing_control_type"  : "step",
            "meshing_frequency"     : 1.0,
            "meshing_before_output" : true,
            "meshing_domains"       : [],
            "write_totalVolumeBeforeMeshing" : true
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.echo_level        = self.settings["echo_level"].GetInt()
        self.dimension         = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]
        self.meshing_frequency = self.settings["meshing_frequency"].GetDouble()
        self.write_total_volume = self.settings["write_totalVolumeBeforeMeshing"].GetBool()

        self.meshing_control_is_time = False
        meshing_control_type   = self.settings["meshing_control_type"].GetString()
        if(meshing_control_type == "time"):
            self.meshing_control_is_time = True
        elif(meshing_control_type == "step"):
            self.meshing_control_is_time = False

        #construct meshing domains
        self.meshing_domains = []
        domains_list = self.settings["meshing_domains"]
        self.number_of_domains = domains_list.size()
        for i in range(0,self.number_of_domains):
            item = domains_list[i]
            domain_module = __import__(item["python_module"].GetString())
            domain = domain_module.CreateMeshingDomain(self.main_model_part,item)
            self.meshing_domains.append(domain)

        # mesh mesher initial values
        self.remesh_domains_active = False
        for domain in self.meshing_domains:
            if( domain.Active() ):
                self.remesh_domains_active = True

        self.neighbours_search_performed = False
        self.step_count   = 1
        self.counter      = 1
        self.next_meshing = 0.0
        self.meshing_before_output = self.settings["meshing_before_output"].GetBool()

    #
    def ExecuteInitialize(self):

        self.fileTotalVolume = None
        self.probe1isolated = None
        self.probe1 = None
        self.probe2 = None
        self.probe3 = None
        self.probe4 = None
        self.probe5 = None
        self.probe6 = None
        self.probe7 = None
        self.probe8 = None
        self.probe9 = None
        self.probe1v = None
        self.probe2v= None
        self.probe3v = None
        self.probe4v = None
        self.probe5v = None
        self.probe6v = None
        self.probe7v = None
        self.probe8v = None
        self.probe9v = None

        print("::[FLUID Meshing_Process]:: meshing frequency", self.meshing_frequency)

        # check restart
        self.restart = False
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.restart = True
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

            if self.meshing_control_is_time:
                self.next_meshing  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] + self.meshing_frequency
            else:
                self.next_meshing = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + self.meshing_frequency
        else:
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
            if self.meshing_control_is_time:
                self.next_meshing  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] + self.meshing_frequency
            else:
                self.next_meshing = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + self.meshing_frequency
            #self.meshing_output = self.meshing_frequency


        self.main_model_part.ProcessInfo.SetValue(KratosDelaunay.INITIALIZED_DOMAINS, False);

        # initialize mesher
        if( self.remesh_domains_active ):

            self.InitializeDomains()

            for domain in self.meshing_domains:
                domain.SetEchoLevel(self.echo_level)
                domain.Initialize()
                if(domain.Active()):
                    domain.ComputeInitialAverageMeshParameters()
                #domain.Check()


    #
    def InitializeDomains(self):

        # initialize the mesher
        if(self.echo_level>1):
            print("::[Remesh_Fluid_Domains]:: Initialize Domains ")

        import domain_utilities
        domain_utils = domain_utilities.DomainUtilities()

        # find node neighbours
        domain_utils.SearchNodeNeighbours(self.main_model_part, self.echo_level)

        # find element neighbours
        domain_utils.SearchElementNeighbours(self.main_model_part, self.echo_level)

        # set neighbour search performed
        self.neighbour_search_performed = True

        # set mesher utilities
        self.mesher_utils = KratosDelaunay.MesherUtilities()

        # set the domain labels to conditions
        self.mesher_utils.SetModelPartNameToConditions(self.main_model_part)

        # find skin and boundary normals
        if(self.restart == False):
            self.BuildMeshBoundaryForFluids()
            #domain_utils.ConstructModelPartBoundary(self.main_model_part, self.echo_level)

            # search nodal h
            if(self.neighbour_search_performed):
                domain_utils.SearchNodalH(self.main_model_part, self.echo_level)

        # set the domain labels to nodes
        self.mesher_utils.SetModelPartNameToNodes(self.main_model_part)

        self.main_model_part.ProcessInfo.SetValue(KratosDelaunay.INITIALIZED_DOMAINS, True)

        if(self.echo_level>1):
            print(self.main_model_part)

    def BuildMeshBoundaryForFluids(self):

        # set building options:

        # define building utility
        model_part_name = self.settings["model_part_name"].GetString()
        #skin_build = KratosDelaunay.BuildModelPartBoundary(self.main_model_part, model_part_name, self.echo_level)
        skin_build = KratosPfemFluid.BuildModelPartBoundaryForFluids(self.main_model_part, model_part_name, self.echo_level)

        # execute building:
        skin_build.Execute()

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1
        currentTime=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        currentStep=self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        if currentStep >= 2 and self.fileTotalVolume is None:
            self.fileTotalVolume = open("totalVolumeBeforeMeshing.txt",'w')
            self.probe1isolated = open("probe1isolated.txt",'w')
            self.probe1 = open("probe1.txt",'w')
            self.probe2 = open("probe2.txt",'w')
            self.probe3 = open("probe3.txt",'w')
            self.probe4 = open("probe4.txt",'w')
            self.probe5 = open("probe5.txt",'w')
            self.probe6 = open("probe6.txt",'w')
            self.probe7 = open("probe7.txt",'w')
            self.probe8 = open("probe8.txt",'w')
            self.probe9 = open("probe9.txt",'w')
            self.probe1v = open("probe1v.txt",'w')
            self.probe2v = open("probe2v.txt",'w')
            self.probe3v = open("probe3v.txt",'w')
            self.probe4v = open("probe4v.txt",'w')
            self.probe5v = open("probe5v.txt",'w')
            self.probe6v = open("probe6v.txt",'w')
            self.probe7v = open("probe7v.txt",'w')
            self.probe8v = open("probe8v.txt",'w')
            self.probe9v = open("probe9v.txt",'w')

        if(currentStep > 1 and self.fileTotalVolume is not None):
            maxYprobe1isolated=0.1
            tolerance=0.015



        if(currentStep > 1 and self.fileTotalVolume is not None):
            maxYprobe1isolated=0.1
            #x1=-0.60
            #x2=-0.65
            #x3=-0.7
            #x4=-0.75
            #x5=-0.8
            #x6=-0.85
            #x7=-0.9
            #x8=-0.95
            #x9=-1.0
            #x1=-1.05
            #x2=-1.1
            #x3=-1.15
            #x4=-1.2
            #x5=-1.25
            #x6=-1.3
            #x7=-1.35
            #x8=-1.4
            #x9=-1.45
            x1=-1.50
            x2=-1.55
            x3=-1.60
            x4=-1.65
            x5=-1.70
            x6=-1.75
            x7=-1.8
            x8=-1.85
            x9=-1.9
            maxYprobe1=-x1*0.577350269
            maxYprobe2=-x2*0.577350269
            maxYprobe3=-x3*0.577350269
            maxYprobe4=-x4*0.577350269
            maxYprobe5=-x5*0.577350269
            maxYprobe6=-x6*0.577350269
            maxYprobe7=-x7*0.577350269
            maxYprobe8=-x8*0.577350269
            maxYprobe9=-x9*0.577350269
            v1=0
            v2=0
            v3=0
            v4=0
            v5=0
            v6=0
            v7=0
            v8=0
            v9=0

            for node in self.main_model_part.Nodes:
                if(node.X>1.87 and node.X<1.93):
                    if(node.Y>maxYprobe1isolated):
                        maxYprobe1isolated=node.Y
                if(node.IsNot(KratosMultiphysics.ISOLATED)):
                    if(node.X>1.87 and node.X<1.93):
                        if(node.Y>maxYprobe1):
                            maxYprobe1=node.Y
                            v1=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>3.07 and node.X<3.13):
                        if(node.Y>maxYprobe2):
                            maxYprobe2=node.Y
                            v2=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>5.97 and node.X<6.03):
                        if(node.Y>maxYprobe3):
                            maxYprobe3=node.Y
                            v3=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>8.77 and node.X<8.83):
                        if(node.Y>maxYprobe4):
                            maxYprobe4=node.Y
                            v4=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>11.27 and node.X<11.33):
                        if(node.Y>maxYprobe5):
                            maxYprobe5=node.Y
                            v5=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>16.77 and node.X<16.83):
                        if(node.Y>maxYprobe6):
                            maxYprobe6=node.Y
                            v6=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>21.97 and node.X<22.03):
                        if(node.Y>maxYprobe7):
                            maxYprobe7=node.Y
                            v7=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>26.57 and node.X<26.63):
                        if(node.Y>maxYprobe8):
                            maxYprobe8=node.Y
                            v8=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)
                    if(node.X>32.47 and node.X<32.53):
                        if(node.Y>maxYprobe9):
                            maxYprobe9=node.Y
                            v9=node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)+node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_y)

           

            outstring = str(currentTime) + " " +  str(maxYprobe1isolated) + "\n"
            self.probe1isolated.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe1) + "\n"
            self.probe1.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe2) + "\n"
            self.probe2.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe3) + "\n"
            self.probe3.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe4) + "\n"
            self.probe4.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe5) + "\n"
            self.probe5.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe6) + "\n"
            self.probe6.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe7) + "\n"
            self.probe7.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe8) + "\n"
            self.probe8.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe9) + "\n"
            self.probe9.write(outstring)
            outstring = str(currentTime) + " " +  str(v1) + "\n"
            self.probe1v.write(outstring)
            outstring = str(currentTime) + " " +  str(v2) + "\n"
            self.probe2v.write(outstring)
            outstring = str(currentTime) + " " +  str(v3) + "\n"
            self.probe3v.write(outstring)
            outstring = str(currentTime) + " " +  str(v4) + "\n"
            self.probe4v.write(outstring)
            outstring = str(currentTime) + " " +  str(v5) + "\n"
            self.probe5v.write(outstring)
            outstring = str(currentTime) + " " +  str(v6) + "\n"
            self.probe6v.write(outstring)
            outstring = str(currentTime) + " " +  str(v7) + "\n"
            self.probe7v.write(outstring)
            outstring = str(currentTime) + " " +  str(v8) + "\n"
            self.probe8v.write(outstring)
            outstring = str(currentTime) + " " +  str(v9) + "\n"
            self.probe9v.write(outstring)
            
            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeAverageMeshParameters()
                    meanVolumeBeforeMeshing=domain.GetMeanVolume()
                    totalVolumeBeforeMeshing=domain.GetTotalVolume()
                    outstring = str(currentTime) + " " +  str(totalVolumeBeforeMeshing) + " "
                    self.fileTotalVolume.write(outstring)
                    #fileTotalVolume = open("totalVolumeBeforeMeshing.txt", 'a')
                    #if(currentStep==2):
                        #fileTotalVolume.seek(0)
                        #fileTotalVolume.truncate()

                    #fileTotalVolume.write(outstring)
                    #fileTotalVolume.close

        volume_acceleration=self.main_model_part.ProcessInfo[KratosMultiphysics.GRAVITY]
        if(currentStep == 1):
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION,volume_acceleration)
  
        if(self.remesh_domains_active or currentStep == 1):
            if( self.meshing_before_output or currentStep == 1):
                if(self.IsMeshingStep() or currentStep == 1):
                    if(self.echo_level>1):
                        print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    print("--> Remesh Fluid Domain at", currentTime, 's')
                    self.RemeshFluidDomains()

        if(currentStep > 1 and self.fileTotalVolume is not None):
            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeAverageMeshParameters()
                    meanVolumeAfterMeshing=domain.GetMeanVolume()
                    totalVolumeAfterMeshing=domain.GetTotalVolume()
                    diffMeanVolume=meanVolumeAfterMeshing-meanVolumeBeforeMeshing
                    diffTotalVolume=totalVolumeAfterMeshing-totalVolumeBeforeMeshing
                    #fileTotalVolume = open("totalVolumeBeforeMeshing.txt", 'a')

                    outstring =  str(totalVolumeAfterMeshing) + " " +  str(diffTotalVolume) + "\n"
                    #fileTotalVolume.write(outstring)
                    #fileTotalVolume.close
                    self.fileTotalVolume.write(outstring)
        if self.fileTotalVolume is not None:
            self.fileTotalVolume.flush()
            self.probe1isolated.flush()
            self.probe1.flush()
            self.probe2.flush()
            self.probe3.flush()
            self.probe4.flush()
            self.probe5.flush()
            self.probe6.flush()
            self.probe7.flush()
            self.probe8.flush()
            self.probe9.flush()
            self.probe1v.flush()
            self.probe2v.flush()
            self.probe3v.flush()
            self.probe4v.flush()
            self.probe5v.flush()
            self.probe6v.flush()
            self.probe7v.flush()
            self.probe8v.flush()
            self.probe9v.flush()

    def ExecuteFinalize(self):
        if self.fileTotalVolume is not None:
            self.fileTotalVolume.close()
            self.probe1isolated.close()
            self.probe1.close()
            self.probe2.close()
            self.probe3.close()
            self.probe4.close()
            self.probe5.close()
            self.probe6.close()
            self.probe7.close()
            self.probe8.close()
            self.probe9.close()
            self.probe1v.close()
            self.probe2v.close()
            self.probe3v.close()
            self.probe4v.close()
            self.probe5v.close()
            self.probe6v.close()
            self.probe7v.close()
            self.probe8v.close()
            self.probe9v.close()
      #if(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] == 1):
          #  for node in self.main_model_part.Nodes:
           #     if (node.Is(KratosMultiphysics.FLUID)):
            #        if(node.IsNot(KratosMultiphysics.RIGID)):
             #           volume_acceleration=node.GetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION)
             #           break
          #  for node in self.main_model_part.Nodes:
            #    if (node.Is(KratosMultiphysics.RIGID)):
             #       node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION,volume_acceleration)


    #
    def ExecuteBeforeOutputStep(self):

        pass
        #if(self.remesh_domains_active):
             #if( self.meshing_before_output ):
                # if(self.IsMeshingStep()):
                   #  print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    # self.RemeshFluidDomains()


    def NodalChecksAndAssignations(self):

        numFluid=0
        numRigid=0
        numRigidFluid=0
        numRigidNotFluid=0
        numBoundary=0
        numIsolated=0
        numFreeSurface=0
        numBlocked=0
        mean_nodal_h=0
        for node in self.main_model_part.Nodes:
            # adding dofs
            if (node.Is(KratosMultiphysics.FLUID)):
                numFluid+=1

                nodal_h=node.GetSolutionStepValue(KratosMultiphysics.NODAL_H)
                mean_nodal_h+=nodal_h

            if (node.Is(KratosMultiphysics.RIGID)):
                numRigid+=1
                if (node.Is(KratosMultiphysics.FLUID)):
                    numRigidFluid+=1
                else:
                    numRigidNotFluid+=1
                    node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0.0)
                if (node.Is(KratosMultiphysics.BOUNDARY)):
                    numBoundary+=1
                if (node.Is(KratosMultiphysics.ISOLATED)):
                    numIsolated+=1
                    node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0.0)

                if (node.Is(KratosMultiphysics.FREE_SURFACE)):
                    numFreeSurface+=1
                if (node.Is(KratosMultiphysics.BLOCKED)):
                    numBlocked+=1

        mean_nodal_h*=1.0/numFluid;
        if(self.echo_level>1):
            print("nodal_h is  ",nodal_h)
            print("numFluid ",numFluid)
            print("numRigid ",numRigid)
            print("numRigidFluid ",numRigidFluid)
            print("numRigidNotFluid ",numRigidNotFluid)
            print("numBoundary ",numBoundary)
            print("numIsolated ",numIsolated)
            print("numFreeSurface ",numFreeSurface)
            print("numBlocked ",numBlocked)

    #
    def ExecuteAfterOutputStep(self):

        if(self.remesh_domains_active):
            if( not self.meshing_before_output ):
                if(self.IsMeshingStep()):
                    self.RemeshFluidDomains()

    #
    def ExecuteMeshing(domain):
        domain.ExecuteMeshing()

    #
    def GetMeshingStep(self):
        return self.counter

    #
    def IsMeshingStep(self):

        if(self.meshing_control_is_time):
            #print( str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME])+">"+ str(self.next_meshing) )
            return ( self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] >= self.next_meshing )
        else:
            return ( self.step_count >= self.next_meshing )

    #

    def RemeshFluidDomains(self):

        if(self.remesh_domains_active):
           # if(self.contact_search):
            #    self.ContactTransfer()

            if( self.echo_level > 1 ):
                print("::[Remesh_fluid_domains_process]:: MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.mesher_utils = KratosDelaunay.MesherUtilities()


            meshing_options.Set(self.mesher_utils.KEEP_ISOLATED_NODES, True)

            #self.model_meshing =  KratosDelaunay.ModelMeshing(self.main_model_part, meshing_options, self.echo_level)
            self.model_meshing =  KratosPfemFluid.ModelMeshingForFluids(self.main_model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()

            id = 0

            for domain in self.meshing_domains:

                domain.ExecuteMeshing();

                self.remesh_executed = True

                id+=1


            self.model_meshing.ExecuteFinalize()

            self.counter += 1
