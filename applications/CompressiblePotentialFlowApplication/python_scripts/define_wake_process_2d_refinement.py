import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication
import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess(Model, settings["Parameters"])


class DefineWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "model_part_name"           : "please specify the model part that contains the kutta nodes",
                "upper_surface_model_part_name" : "please specify the model part that contains the upper surface nodes",
                "lower_surface_model_part_name" : "please specify the model part that contains the lower surface nodes",
                "fluid_part_name"           : "MainModelPart",
                "direction"                 : [1.0,0.0,0.0],
                "stl_filename"              : "please specify name of stl file",
                "epsilon"    : 1e-9
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings)

        self.direction = KratosMultiphysics.Vector(3)
        self.direction[0] = settings["direction"][0].GetDouble()
        self.direction[1] = settings["direction"][1].GetDouble()
        self.direction[2] = settings["direction"][2].GetDouble()
        dnorm = math.sqrt(self.direction[0]**2 + self.direction[1]**2 + self.direction[2]**2)
        self.direction[0] /= dnorm
        self.direction[1] /= dnorm
        self.direction[2] /= dnorm
        print(self.direction)

        self.epsilon = settings["epsilon"].GetDouble()

        self.upper_surface_model_part = Model[settings["upper_surface_model_part_name"].GetString()]
        self.lower_surface_model_part = Model[settings["lower_surface_model_part_name"].GetString()]
        self.fluid_model_part = Model[settings["fluid_part_name"].GetString()]
        self.wake_model_part = self.fluid_model_part.CreateSubModelPart("wake_modelpart")
        self.kutta_model_part = self.fluid_model_part.CreateSubModelPart("kutta_model_part")
        self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_model_part")
        self.penalty_model_part = self.fluid_model_part.CreateSubModelPart("penalty_model_part")

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part,
                                                                       self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_model_part, AvgElemNum, AvgNodeNum)
        # Find neighbours
        nodal_neighbour_search.Execute()

        self.stl_filename = settings["stl_filename"].GetString()
        
    def Execute(self):
        #find the trailing edge node
        pos = -1000
        pos_up = -1000
        pos_down = -1000
        
        for unode in self.upper_surface_model_part.Nodes:
            if(unode.X > pos):
                pos = unode.X
                te_node = unode

        for unode in self.upper_surface_model_part.Nodes:
            if(unode.X > pos_up and unode.X < pos - self.epsilon):
                pos_up = unode.X
                upper_te_node = unode

        self.kutta_model_part.Nodes.append(te_node)
        print('kutta node = ', te_node)
        print('upper trailing edge node =', upper_te_node)
        
        for lnode in self.lower_surface_model_part.Nodes:
            if(lnode.X > pos_down and lnode.X < pos -self.epsilon):
                pos_down = lnode.X
                lower_te_node = lnode

        print('lower trailing edge node =', lower_te_node)
        print('...Finding the trailing edge finished...')

        #Create trailing edge element
        prop = self.fluid_model_part.GetProperties()[1]
        elem_id = self.fluid_model_part.NumberOfElements() + 1
        self.fluid_model_part.CreateNewElement("CompressiblePotentialFlowElement2D3N", elem_id, [te_node.Id, upper_te_node.Id, lower_te_node.Id], prop)

        te_elem = self.fluid_model_part.GetElements()[elem_id]
        te_elem.Set(KratosMultiphysics.MODIFIED, True)
        te_elem.Set(KratosMultiphysics.STRUCTURE, True)

        #mark as STRUCTURE and deactivate the elements that touch the kutta node
        for node in self.kutta_model_part.Nodes:
            node.Set(KratosMultiphysics.STRUCTURE)
            node.Set(KratosMultiphysics.THERMAL)
            x1 = node.X
            y1 = node.Y
            z1 = node.Z

        #Selecting upper and lower surface nodes
        for node in self.upper_surface_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.UPPER_SURFACE, True)
        for node in self.lower_surface_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.LOWER_SURFACE, True)
        
        ##find wake node in the outflow boundary
        #pos = 0
        #print(self.kutta_model_part)

        #compute the distances of the elements of the wake, and decide which ones are wake
        xn = KratosMultiphysics.Vector(3)

        self.n = KratosMultiphysics.Vector(3)
        self.n[0] = -self.direction[1]
        self.n[1] = self.direction[0]
        self.n[2] = 0.0

        for node in self.kutta_model_part.Nodes:
            x0 = node.X
            y0 = node.Y
            for elem in self.fluid_model_part.Elements:
                #check in the potentially active portion
                potentially_active_portion = False
                for elnode in elem.GetNodes():
                    #all nodes that touch the kutta nodes are potentiallyactive
                    if(elnode.Is(KratosMultiphysics.STRUCTURE)):
                        potentially_active_portion = True
                        self.trailing_edge_model_part.Elements.append(elem)
                        self.penalty_model_part.Elements.append(elem)
                        #elem.Set(KratosMultiphysics.THERMAL)
                        break
                for elnode in elem.GetNodes():
                    #compute x distance
                    xn[0] = elnode.X - x0
                    xn[1] = elnode.Y - y0
                    xn[2] = 0.0
                    dx = xn[0]*self.direction[0] + xn[1]*self.direction[1]
                    if(dx > 0):
                        potentially_active_portion = True
                        break

                if(potentially_active_portion):
                    distances = KratosMultiphysics.Vector(len(elem.GetNodes()))

                    counter = 0
                    for elnode in elem.GetNodes():
                        xn[0] = elnode.X - x0
                        xn[1] = elnode.Y - y0
                        xn[2] = 0.0
                        d = xn[0]*self.n[0] + xn[1]*self.n[1]
                        if(abs(d) < self.epsilon):
                            d = self.epsilon
                        if(d<0 and
                            elnode.IsNot(KratosMultiphysics.STRUCTURE) and
                            elnode.GetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.UPPER_SURFACE) == True):
                            d = self.epsilon
                            print('\ndetected upper surface node')
                            print(elnode)
                        if( xn[0] < 0 and 
                            d > 0 and #for high angles of attack (selecting nodes in the lower surface)
                            elnode.IsNot(KratosMultiphysics.STRUCTURE) and
                            elnode.GetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.LOWER_SURFACE) == True):
                            d = -self.epsilon
                            print('\ndetected lower surface node!')
                            print(elnode)
                        distances[counter] = d
                        counter += 1

                    npos = 0
                    nneg = 0
                    for d in distances:
                        if(d < 0):
                            nneg += 1
                        else:
                            npos += 1

                    if(nneg > 0 and npos > 0):
                        elem.Set(KratosMultiphysics.MARKER, True)
                        self.wake_model_part.Elements.append(elem)
                        counter = 0
                        for elnode in elem.GetNodes():
                            elnode.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distances[counter])
                            elnode.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,5.0)
                            counter += 1
                            #if(elnode.Is(KratosMultiphysics.STRUCTURE)):
                            #    #selecting Kutta elements
                            #    elem.Set(KratosMultiphysics.STRUCTURE)
                            #    elem.Set(KratosMultiphysics.MODIFIED)
                            '''
                            dx = elnode.X - x1
                            dy = elnode.Y - y1
                            dz = elnode.Z - z1
                            tmp = dx * \
                                self.direction[0] + dy * \
                                self.direction[1] + dz*self.direction[2]
                            if(tmp > pos):
                               pos = tmp
                            '''
                        elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, distances)

                        #for elnode in elem.GetNodes():
                        #if elnode.Is(KratosMultiphysics.STRUCTURE):
                        #elem.Set(KratosMultiphysics.ACTIVE,False)
                        #elem.Set(KratosMultiphysics.MARKER,False)
        
        #'''
        #mark only airfoil elements as kutta elements
        for trailing_elem in self.trailing_edge_model_part.Elements:
            up_counter = 0
            low_counter = 0
            for wake_nodes in trailing_elem.GetNodes():
                for upnode in self.upper_surface_model_part.Nodes:
                    if(abs(wake_nodes.X-upnode.X) < 1e-6 and abs(wake_nodes.Y-upnode.Y) < 1e-6):
                        up_counter += 1
                for lownode in self.lower_surface_model_part.Nodes:
                    if(abs(wake_nodes.X-lownode.X) < 1e-6 and abs(wake_nodes.Y-lownode.Y) < 1e-6):
                        low_counter += 1
            if(low_counter > 1):# or low_counter > 1):
                #trailing_elem.Set(KratosMultiphysics.STRUCTURE)
                #trailing_elem.Set(KratosMultiphysics.MODIFIED, True)
                print('modified element =', trailing_elem)
                for node in trailing_elem.GetNodes():
                    node.Set(KratosMultiphysics.STRUCTURE)
        #'''

        # '''
        #loop to select as many modified elements as desired
        extra_penalty_elements_layers = 2
        for i in range(extra_penalty_elements_layers):
            for elem in self.penalty_model_part.Elements:
                for elnode in elem.GetNodes():
                    elnode.Set(KratosMultiphysics.THERMAL)

            for elem in self.fluid_model_part.Elements:
                counter = 0
                for elnode in elem.GetNodes():
                    if(elnode.Is(KratosMultiphysics.THERMAL)):
                        counter += 1
                if(counter > 0):
                    elem.Set(KratosMultiphysics.THERMAL)
                    self.penalty_model_part.Elements.append(elem)
        #'''

        #'''
        #tolerance = 1e-5
        counter_non_marked_elements = 0
        counter_marked_elements = 0
        for elem in self.trailing_edge_model_part.Elements:
            if(elem.Is(KratosMultiphysics.MARKER)):
                counter_marked_elements +=1
            if(elem.IsNot(KratosMultiphysics.MARKER)):
                counter_non_marked_elements += 1
                '''
                #elem.Set(KratosMultiphysics.MARKER, True)
                #elem.Set(KratosMultiphysics.ISOLATED, True)
                distances = KratosMultiphysics.Vector(len(elem.GetNodes()))
                counter = 0
                for elnode in elem.GetNodes():
                    xn[0] = elnode.X - x0
                    xn[1] = elnode.Y - y0
                    xn[2] = 0.0
                    d = xn[0]*self.n[0] + xn[1]*self.n[1]
                    #if( abs(d) < tolerance and elnode.IsNot(KratosMultiphysics.STRUCTURE)):
                    #    d = -d
                    if( elnode.Is(KratosMultiphysics.STRUCTURE)):
                        d = self.epsilon
                    distances[counter] = d
                    counter += 1
                npos = 0
                nneg = 0
                for d in distances:
                    if(d < 0):
                        nneg += 1
                    else:
                        npos += 1

                if(nneg > 0 and npos > 0):
                    print(elem)
                    print(distances)
                    #elem.Set(KratosMultiphysics.MARKER, True)
                    counter = 0
                    for elnode in elem.GetNodes():
                        elnode.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distances[counter])
                        counter += 1
                elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, distances)
        #'''
        '''
        if(counter_non_marked_elements > 2 or counter_marked_elements > 3):
            print('counter_non_marked_elements = ', counter_non_marked_elements)
            print('counter_marked_elements = ', counter_marked_elements)
            extra_kutta_elements_number = 1
        else:
            extra_kutta_elements_number = 0
        '''
        
        '''
        #loop to select as many kutta elements as desired
        extra_kutta_elements_number = 2
        for i in range(extra_kutta_elements_number):
            for elem in self.kutta_model_part.Elements:
                if(elem.Is(KratosMultiphysics.MARKER)):
                    for elnode in elem.GetNodes():
                        elnode.Set(KratosMultiphysics.STRUCTURE)

            for elem in self.wake_model_part.Elements:
                counter = 0
                for elnode in elem.GetNodes():
                    if(elnode.Is(KratosMultiphysics.STRUCTURE)):
                        counter += 1
                if(counter > 1):
                    elem.Set(KratosMultiphysics.STRUCTURE)
                    self.kutta_model_part.Elements.append(elem)
        #'''

        '''
        #find the wake elements at the outflow boundary
        for elem in self.wake_model_part.Elements:
            counter = 0
            for elnode in elem.GetNodes():
                dx = elnode.X - x1
                dy = elnode.Y - y1
                dz = elnode.Z - z1
                tmp = dx*self.direction[0] + dy*self.direction[1] + dz*self.direction[2]
                
                if(tmp > pos/4e4-1e-9):
                    counter +=1
            if(counter > 1):
                elem.Set(KratosMultiphysics.BOUNDARY,True)
                for elnode in elem.GetNodes():
                    elnode.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,10.0)
                    #node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,d[counter])
                #print('outlet boundary element found')
                #print(elem)
                #print(tmp)
                    
        #'''

    def ExecuteInitialize(self):
        self.Execute()
