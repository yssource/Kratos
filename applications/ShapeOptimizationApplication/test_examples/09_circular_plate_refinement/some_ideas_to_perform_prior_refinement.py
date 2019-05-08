
# --------------------------------------------------------------------------
    def PerformAPrioriRefinement(self):

        print("\n> Initializing prior refinement...")
        start_time = time.time()

        max_angular_change = 2 # in degree
        non_linearity_threshold = 0.2

        temp_cf = ConditionsFactory(self.fe_model_part, self.cad_model, self.parameters)
        boundary_polygons = temp_cf.GetBoundaryPolygons()
        boundary_tessellation_tolerance = temp_cf.boundary_tessellation_tolerance

        # Collect points to be mapped
        temp_model = KratosMultiphysics.Model()
        destination_mdpa = temp_model.CreateModelPart("temp_model_part")
        destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE_ABSOLUTE)
        destination_mdpa.AddNodalSolutionStepVariable(KratosShape.GRAD_SHAPE_CHANGE)

        node_itr = 0

        list_of_parameter_positions_inside = {face.Key(): [] for face in self.cad_model.GetByType('BrepFace')}
        id_to_itr = {}

        # Prepare mapping
        for face in self.cad_model.GetByType('BrepFace'):
            geometry = face.Data().Geometry()
            geometry_data = geometry.Data()

            # Identify elements to refine
            intervals_along_u_to_refine = [[]]

            knots_u = geometry_data.KnotsU()
            knots_v = geometry_data.KnotsV()

            # Remove duplicated values
            u_list = list(set(knots_u))
            v_list = list(set(knots_v))
            u_list.sort()
            v_list.sort()

            for v_itr in range(len(v_list)-1):
                for u_itr in range(len(u_list)-1):
                    u = (u_list[u_itr+1] + u_list[u_itr]) / 2
                    v = (v_list[v_itr+1] + v_list[v_itr]) / 2

                    is_inside, is_on_boundary = temp_cf.Contains((u, v), boundary_polygons[face.Key()], boundary_tessellation_tolerance*1.1)
                    if is_inside:

                        point = geometry_data.PointAt(u, v)
                        # point_ptr = self.cad_model.Add(an.Point3D(location=point))
                        # point_ptr.Attributes().SetLayer("CenterPoints_"+str(face.Key()))

                        destination_mdpa.CreateNewNode(node_itr+1, point[0], point[1], point[2])

                        id_to_itr[(face.Key(), u, v)] = node_itr

                        list_of_parameter_positions_inside[face.Key()].append((u,v))

                        node_itr += 1

        # Map information from fem to integration points using element based mapper
        mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, self.parameters["conditions"]["general"]["mapping_cad_fem"].Clone() )
        mapper.Map( KratosShape.GRAD_SHAPE_CHANGE, KratosShape.GRAD_SHAPE_CHANGE )

        # Extract info from new fe model part
        list_of_nodes = [node for node in destination_mdpa.Nodes]
        list_of_grads = [np.array(node.GetSolutionStepValue(KratosShape.GRAD_SHAPE_CHANGE)) for node in destination_mdpa.Nodes]

        l2_norms_grads = np.array([ np.linalg.norm(grad) for grad in list_of_grads])
        max_norm_grads = max(l2_norms_grads)

        # Identify spots to refine
        intervals_along_u_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}
        intervals_along_v_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}

        for face in self.cad_model.GetByType('BrepFace'):
            geometry = face.Data().Geometry()
            geometry_data = geometry.Data()

            relevant_parameters = list_of_parameter_positions_inside[face.Key()]

            for (u,v) in relevant_parameters:

                        node = list_of_nodes[id_to_itr[(face.Key(),u,v)]]
                        gradient = list_of_grads[id_to_itr[(face.Key(),u,v)]]

                        if np.linalg.norm(gradient) > 1:
                            pv = geometry_data.PointAt(u, v)
                            point_ptr = self.cad_model.Add(an.Point3D(location=pv))
                            point_ptr.Attributes().SetLayer("PositionsToRefine")

        self.__SaveCadModel("cad_after_prior_refinement.iga")

        print("> Prior refinement finished in" ,round( time.time()-start_time, 3 ), " s.")

        return


# --------------------------------------------------------------------------
    def PerformAPrioriRefinement(self):

        print("\n> Initializing prior refinement...")
        start_time = time.time()

        max_angular_change = 2 # in degree
        non_linearity_threshold = 0.2

        temp_cf = ConditionsFactory(self.fe_model_part, self.cad_model, self.parameters)
        boundary_polygons = temp_cf.GetBoundaryPolygons()
        boundary_tessellation_tolerance = temp_cf.boundary_tessellation_tolerance

        # Collect points to be mapped
        temp_model = KratosMultiphysics.Model()
        destination_mdpa = temp_model.CreateModelPart("temp_model_part")
        destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE_ABSOLUTE)
        destination_mdpa.AddNodalSolutionStepVariable(KratosShape.GRAD_SHAPE_CHANGE)

        node_itr = 0

        list_of_parameter_positions_inside = {face.Key(): [] for face in self.cad_model.GetByType('BrepFace')}
        id_to_itr = {}

        # Prepare mapping
        for face in self.cad_model.GetByType('BrepFace'):
            geometry = face.Data().Geometry()
            geometry_data = geometry.Data()

            # Identify elements to refine
            intervals_along_u_to_refine = [[]]

            knots_u = geometry_data.KnotsU()
            knots_v = geometry_data.KnotsV()

            # Remove duplicated values
            u_list = list(set(knots_u))
            v_list = list(set(knots_v))
            u_list.sort()
            v_list.sort()

            for v in v_list:
                for u in u_list:
                    is_inside, is_on_boundary = temp_cf.Contains((u, v), boundary_polygons[face.Key()], boundary_tessellation_tolerance*1.1)
                    if is_inside:

                        point = geometry_data.PointAt(u, v)
                        # point_ptr = self.cad_model.Add(an.Point3D(location=point))
                        # point_ptr.Attributes().SetLayer("CenterPoints_"+str(face.Key()))

                        destination_mdpa.CreateNewNode(node_itr+1, point[0], point[1], point[2])

                        id_to_itr[(face.Key(), u, v)] = node_itr

                        list_of_parameter_positions_inside[face.Key()].append((u,v))

                        node_itr += 1

        # Map information from fem to integration points using element based mapper
        mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, self.parameters["conditions"]["general"]["mapping_cad_fem"].Clone() )
        mapper.Map( KratosShape.SHAPE_CHANGE_ABSOLUTE, KratosShape.SHAPE_CHANGE_ABSOLUTE )
        mapper.Map( KratosShape.GRAD_SHAPE_CHANGE, KratosShape.GRAD_SHAPE_CHANGE )

        # Extract info from new fe model part
        list_of_nodes = [node for node in destination_mdpa.Nodes]
        list_of_updates = [node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE_ABSOLUTE) for node in destination_mdpa.Nodes]
        list_of_grads = [np.array(node.GetSolutionStepValue(KratosShape.GRAD_SHAPE_CHANGE)) for node in destination_mdpa.Nodes]

        # Identify spots to refine
        intervals_along_u_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}
        intervals_along_v_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}

        l2_norms = np.array([ np.linalg.norm(update) for update in list_of_updates])
        max_norm = max(l2_norms)

        for face in self.cad_model.GetByType('BrepFace'):
            geometry = face.Data().Geometry()
            geometry_data = geometry.Data()

            knots_u = geometry_data.KnotsU()
            knots_v = geometry_data.KnotsV()

            # Remove duplicated values
            u_list = list(set(knots_u))
            v_list = list(set(knots_v))
            u_list.sort()
            v_list.sort()

            relevant_parameters = list_of_parameter_positions_inside[face.Key()]

            for v_itr, v in enumerate(v_list):
                for u_itr, u in enumerate(u_list):
                    if (u,v) in relevant_parameters:

                        node = list_of_nodes[id_to_itr[(face.Key(),u,v)]]

                        point_00 = np.array([node.X, node.Y, node.Z])
                        update_00 = list_of_updates[id_to_itr[(face.Key(),u,v)]]
                        gradient_00 = list_of_grads[id_to_itr[(face.Key(),u,v)]]

                            # point_ptr = self.cad_model.Add(an.Point3D(location=point_00))
                            # point_ptr.Attributes().SetLayer("PointsWithUpdatesBigger_7.2"+str(face.Key()))

                        u1 = u_list[u_itr+1]
                        v1 = v_list[v_itr+1]

                        if (u1,v) in relevant_parameters:
                            point_10 = geometry_data.PointAt(u1, v)
                            direction = point_10-point_00
                            distance = np.linalg.norm(direction)
                            normalized_direction = direction/distance

                            update_10 = list_of_updates[id_to_itr[(face.Key(),u1,v)]]
                            delta_update = update_00-update_10

                            directional_gradient = np.dot(gradient_00, normalized_direction)

                            nonlinearity_measure = (update_00 + directional_gradient*distance - update_10)


                            if abs(nonlinearity_measure) > 1:
                                print("u--")
                                print(nonlinearity_measure)
                                print(update_00)
                                print("target = ", update_10)
                                print("estimate = ", update_00 + directional_gradient*distance)

                                pv = geometry_data.PointAt(u1, v)
                                point_ptr = self.cad_model.Add(an.Point3D(location=pv))
                                point_ptr.Attributes().SetLayer("UsToRefine")


                        if (u,v1) in relevant_parameters:
                            point_01 = geometry_data.PointAt(u, v1)
                            direction = point_01-point_00
                            distance = np.linalg.norm(direction)
                            normalized_direction = direction/distance

                            update_01 = list_of_updates[id_to_itr[(face.Key(),u,v1)]]
                            delta_update = update_00-update_01

                            directional_gradient = np.dot(gradient_00, normalized_direction)

                            nonlinearity_measure = (update_00 + directional_gradient*distance - update_01)


                            # if nonlinearity_measure > non_linearity_threshold:
                            #     exist_intervals_to_refine = self.__AddVIntervalToList(geometry, u, (v1+v)*0.5, intervals_along_u_to_refine)

                            if abs(nonlinearity_measure) > 1:
                                print("v--")
                                print(nonlinearity_measure)
                                print(update_00)
                                print(update_01)

                                pv = geometry_data.PointAt(u, v1)
                                point_ptr = self.cad_model.Add(an.Point3D(location=pv))
                                point_ptr.Attributes().SetLayer("VsToRefine")


                        node_itr += 1


        self.__SaveCadModel("cad_after_prior_refinement.iga")

        print("> Prior refinement finished in" ,round( time.time()-start_time, 3 ), " s.")

        return

    # # --------------------------------------------------------------------------
    # def PerformAPrioriRefinement(self):

    #     print("\n> Initializing prior refinement...")
    #     start_time = time.time()

    #     max_angular_change = 2 # in degree
    #     non_linearity_threshold = 0.2

    #     temp_cf = ConditionsFactory(self.fe_model_part, self.cad_model, self.parameters)
    #     boundary_polygons = temp_cf.GetBoundaryPolygons()
    #     boundary_tessellation_tolerance = temp_cf.boundary_tessellation_tolerance

    #     # Collect points to be mapped
    #     temp_model = KratosMultiphysics.Model()
    #     destination_mdpa = temp_model.CreateModelPart("temp_model_part")
    #     destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

    #     node_itr = 0

    #     list_of_parameter_positions_inside = {face.Key(): [] for face in self.cad_model.GetByType('BrepFace')}
    #     id_to_itr = {}

    #     # Prepare mapping
    #     for face in self.cad_model.GetByType('BrepFace'):
    #         geometry = face.Data().Geometry()
    #         geometry_data = geometry.Data()

    #         # Identify elements to refine
    #         intervals_along_u_to_refine = [[]]

    #         knots_u = geometry_data.KnotsU()
    #         knots_v = geometry_data.KnotsV()

    #         # Remove duplicated values
    #         u_list = list(set(knots_u))
    #         v_list = list(set(knots_v))
    #         u_list.sort()
    #         v_list.sort()

    #         for v in v_list:
    #             for u in u_list:
    #                 is_inside, is_on_boundary = temp_cf.Contains((u, v), boundary_polygons[face.Key()], boundary_tessellation_tolerance*1.1)
    #                 if is_inside:

    #                     point = geometry_data.PointAt(u, v)
    #                     # point_ptr = self.cad_model.Add(an.Point3D(location=point))
    #                     # point_ptr.Attributes().SetLayer("CenterPoints_"+str(face.Key()))

    #                     destination_mdpa.CreateNewNode(node_itr+1, point[0], point[1], point[2])

    #                     id_to_itr[(face.Key(), u, v)] = node_itr

    #                     list_of_parameter_positions_inside[face.Key()].append((u,v))

    #                     node_itr += 1

    #     # Map information from fem to integration points using element based mapper
    #     mapper_parameters = KratosMultiphysics.Parameters("""{
    #         "mapper_type" : "nearest_element",
    #         "search_radius" : 1.0
    #     }""")
    #     mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, mapper_parameters )
    #     mapper.Map( KratosShape.SHAPE_CHANGE, KratosShape.SHAPE_CHANGE )

    #     # Extract info from new fe model part
    #     list_of_nodes = [node for node in destination_mdpa.Nodes]
    #     list_of_updates = [np.array(node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)) for node in destination_mdpa.Nodes]

    #     # Identify spots to refine
    #     intervals_along_u_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}
    #     intervals_along_v_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}

    #     l2_norms = np.array([ np.linalg.norm(update) for update in list_of_updates])
    #     max_norm = max(l2_norms)

    #     for face in self.cad_model.GetByType('BrepFace'):
    #         geometry = face.Data().Geometry()
    #         geometry_data = geometry.Data()

    #         knots_u = geometry_data.KnotsU()
    #         knots_v = geometry_data.KnotsV()

    #         # Remove duplicated values
    #         u_list = list(set(knots_u))
    #         v_list = list(set(knots_v))
    #         u_list.sort()
    #         v_list.sort()

    #         relevant_parameters = list_of_parameter_positions_inside[face.Key()]

    #         for v_itr, v in enumerate(v_list):
    #             for u_itr, u in enumerate(u_list):
    #                 if (u,v) in relevant_parameters:

    #                     node = list_of_nodes[id_to_itr[(face.Key(),u,v)]]

    #                     point_00 = np.array([node.X, node.Y, node.Z])
    #                     update_00 = np.array(node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE))

    #                         # point_ptr = self.cad_model.Add(an.Point3D(location=point_00))
    #                         # point_ptr.Attributes().SetLayer("PointsWithUpdatesBigger_7.2"+str(face.Key()))

    #                     u1 = u_list[u_itr+1]
    #                     v1 = v_list[v_itr+1]

    #                     if (u1,v) in relevant_parameters:
    #                         point_10 = geometry_data.PointAt(u1, v)
    #                         # distance = np.linalg.norm(point_10-point_00)

    #                         update_10 = list_of_updates[id_to_itr[(face.Key(),u1,v)]]
    #                         delta_update = np.linalg.norm(update_00-update_10)

    #                         nonlinearity_measure = delta_update / max_norm

    #                         print(nonlinearity_measure)

    #                         if nonlinearity_measure > non_linearity_threshold:
    #                             exist_intervals_to_refine = self.__AddUIntervalToList(geometry, (u1+u)*0.5, v, intervals_along_u_to_refine)

    #                             pu = geometry_data.PointAt((u1+u)*0.5, v)
    #                             point_ptr = self.cad_model.Add(an.Point3D(location=pu))
    #                             point_ptr.Attributes().SetLayer("UsToRefine")


    #                     if (u,v1) in relevant_parameters:
    #                         point_01 = geometry_data.PointAt(u, v1)
    #                         # distance = np.linalg.norm(point_01-point_00)

    #                         update_01 = list_of_updates[id_to_itr[(face.Key(),u,v1)]]
    #                         delta_update = np.linalg.norm(update_00-update_01)

    #                         nonlinearity_measure = delta_update / max_norm

    #                         print(nonlinearity_measure)

    #                         if nonlinearity_measure > non_linearity_threshold:
    #                             exist_intervals_to_refine = self.__AddVIntervalToList(geometry, u, (v1+v)*0.5, intervals_along_u_to_refine)


    #                             pv = geometry_data.PointAt(u, (v1+v)*0.5)
    #                             point_ptr = self.cad_model.Add(an.Point3D(location=pv))
    #                             point_ptr.Attributes().SetLayer("VsToRefine")


    #                     node_itr += 1


    #     self.__SaveCadModel("cad_after_prior_refinement.iga")

    #     print("> Prior refinement finished in" ,round( time.time()-start_time, 3 ), " s.")


    #     err
    #     return

    # # --------------------------------------------------------------------------
    # def PerformAPrioriRefinement(self):

    #     print("\n> Initializing prior refinement...")
    #     start_time = time.time()

    #     max_angular_change = 2 # in degree
    #     non_linearity_threshold = 0.2

    #     temp_cf = ConditionsFactory(self.fe_model_part, self.cad_model, self.parameters)
    #     boundary_polygons = temp_cf.GetBoundaryPolygons()
    #     boundary_tessellation_tolerance = temp_cf.boundary_tessellation_tolerance

    #     # Collect points to be mapped
    #     temp_model = KratosMultiphysics.Model()
    #     destination_mdpa = temp_model.CreateModelPart("temp_model_part")
    #     destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

    #     node_itr = 0

    #     list_of_parameter_positions_inside = {face.Key(): [] for face in self.cad_model.GetByType('BrepFace')}
    #     id_to_itr = {}

    #     # Prepare mapping
    #     for face in self.cad_model.GetByType('BrepFace'):
    #         geometry = face.Data().Geometry()
    #         geometry_data = geometry.Data()

    #         # Identify elements to refine
    #         intervals_along_u_to_refine = [[]]

    #         knots_u = geometry_data.KnotsU()
    #         knots_v = geometry_data.KnotsV()

    #         # Remove duplicated values
    #         u_list = list(set(knots_u))
    #         v_list = list(set(knots_v))
    #         u_list.sort()
    #         v_list.sort()

    #         for v in v_list:
    #             for u in u_list:
    #                 is_inside, is_on_boundary = temp_cf.Contains((u, v), boundary_polygons[face.Key()], boundary_tessellation_tolerance*1.1)
    #                 if is_inside:

    #                     point = geometry_data.PointAt(u, v)
    #                     # point_ptr = self.cad_model.Add(an.Point3D(location=point))
    #                     # point_ptr.Attributes().SetLayer("CenterPoints_"+str(face.Key()))

    #                     destination_mdpa.CreateNewNode(node_itr+1, point[0], point[1], point[2])

    #                     id_to_itr[(face.Key(), u, v)] = node_itr

    #                     list_of_parameter_positions_inside[face.Key()].append((u,v))

    #                     node_itr += 1

    #     # Map information from fem to integration points using element based mapper
    #     mapper_parameters = KratosMultiphysics.Parameters("""{
    #         "mapper_type" : "nearest_element",
    #         "search_radius" : 1.0
    #     }""")
    #     mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, mapper_parameters )
    #     mapper.Map( KratosShape.SHAPE_CHANGE, KratosShape.SHAPE_CHANGE )

    #     # Extract info from new fe model part
    #     list_of_nodes = [node for node in destination_mdpa.Nodes]
    #     list_of_updates = [np.array(node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)) for node in destination_mdpa.Nodes]

    #     # Identify spots to refine
    #     intervals_along_u_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}
    #     intervals_along_v_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}

    #     l2_norms = np.array([ np.linalg.norm(update) for update in list_of_updates])
    #     max_norm = max(l2_norms)

    #     for face in self.cad_model.GetByType('BrepFace'):
    #         geometry = face.Data().Geometry()
    #         geometry_data = geometry.Data()

    #         knots_u = geometry_data.KnotsU()
    #         knots_v = geometry_data.KnotsV()

    #         # Remove duplicated values
    #         u_list = list(set(knots_u))
    #         v_list = list(set(knots_v))
    #         u_list.sort()
    #         v_list.sort()

    #         relevant_parameters = list_of_parameter_positions_inside[face.Key()]

    #         for v_itr, v in enumerate(v_list):
    #             for u_itr, u in enumerate(u_list):
    #                 if (u,v) in relevant_parameters:

    #                     node = list_of_nodes[id_to_itr[(face.Key(),u,v)]]

    #                     point_00 = np.array([node.X, node.Y, node.Z])
    #                     update_00 = np.array(node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE))

    #                         # point_ptr = self.cad_model.Add(an.Point3D(location=point_00))
    #                         # point_ptr.Attributes().SetLayer("PointsWithUpdatesBigger_7.2"+str(face.Key()))

    #                     u1 = u_list[u_itr+1]
    #                     v1 = v_list[v_itr+1]

    #                     if (u1,v) in relevant_parameters:
    #                         point_10 = geometry_data.PointAt(u1, v)
    #                         # distance = np.linalg.norm(point_10-point_00)

    #                         update_10 = list_of_updates[id_to_itr[(face.Key(),u1,v)]]
    #                         delta_update = np.linalg.norm(update_00-update_10)

    #                         nonlinearity_measure = delta_update / max_norm

    #                         print(nonlinearity_measure)

    #                         if nonlinearity_measure > non_linearity_threshold:
    #                             exist_intervals_to_refine = self.__AddUIntervalToList(geometry, (u1+u)*0.5, v, intervals_along_u_to_refine)

    #                             pu = geometry_data.PointAt((u1+u)*0.5, v)
    #                             point_ptr = self.cad_model.Add(an.Point3D(location=pu))
    #                             point_ptr.Attributes().SetLayer("UsToRefine")


    #                     if (u,v1) in relevant_parameters:
    #                         point_01 = geometry_data.PointAt(u, v1)
    #                         # distance = np.linalg.norm(point_01-point_00)

    #                         update_01 = list_of_updates[id_to_itr[(face.Key(),u,v1)]]
    #                         delta_update = np.linalg.norm(update_00-update_01)

    #                         nonlinearity_measure = delta_update / max_norm

    #                         print(nonlinearity_measure)

    #                         if nonlinearity_measure > non_linearity_threshold:
    #                             exist_intervals_to_refine = self.__AddVIntervalToList(geometry, u, (v1+v)*0.5, intervals_along_u_to_refine)


    #                             pv = geometry_data.PointAt(u, (v1+v)*0.5)
    #                             point_ptr = self.cad_model.Add(an.Point3D(location=pv))
    #                             point_ptr.Attributes().SetLayer("VsToRefine")


    #                     node_itr += 1


    #     self.__SaveCadModel("cad_after_prior_refinement.iga")

    #     print("> Prior refinement finished in" ,round( time.time()-start_time, 3 ), " s.")


    #     err
    #     return



    # First draft
    # # --------------------------------------------------------------------------
    # def PerformAPrioriRefinement(self):

    #     print("\n> Initializing prior refinement...")
    #     start_time = time.time()


    #     refinement_tolerance = 0.1

    #     temp_cf = ConditionsFactory(self.fe_model_part, self.cad_model, self.parameters)
    #     boundary_polygons = temp_cf.GetBoundaryPolygons()
    #     boundary_tessellation_tolerance = temp_cf.boundary_tessellation_tolerance

    #     # Collect points to be mapped
    #     temp_model = KratosMultiphysics.Model()
    #     destination_mdpa = temp_model.CreateModelPart("temp_model_part")
    #     destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

    #     node_itr = 1

    #     id_to_itr = {}
    #     relevant_us_for_given_vs = {face.Key(): []} for face in self.cad_model.GetByType('BrepFace')}

    #     # Perform actual refinement
    #     for face in self.cad_model.GetByType('BrepFace'):
    #         geometry = face.Data().Geometry()
    #         geometry_data = geometry.Data()

    #         # Identify elements to refine
    #         intervals_along_u_to_refine = [[]]

    #         knots_u = geometry_data.KnotsU()
    #         knots_v = geometry_data.KnotsV()

    #         # Remove duplicated values
    #         u_list = list(set(knots_u))
    #         v_list = list(set(knots_v))
    #         u_list.sort()
    #         v_list.sort()

    #         for v in v_list:
    #             relevant_us_for_given_vs[face.Key()].append([])
    #             for u in u_list:
    #                 is_inside, is_on_boundary = temp_cf.Contains((u, v), boundary_polygons[face.Key()], boundary_tessellation_tolerance*1.1)
    #                 if is_inside:

    #                     point = geometry_data.PointAt(u, v)
    #                     point_ptr = self.cad_model.Add(an.Point3D(location=point))
    #                     point_ptr.Attributes().SetLayer("CenterPoints_"+str(face.Key()))

    #                     destination_mdpa.CreateNewNode(node_itr, point[0], point[1], point[2])

    #                     relevant_us_for_given_vs[face.Key()[:].append(u)
    #                     id_to_itr[(face.Key(), u, v)] = node_itr

    #                     node_itr += 1

    #     # Map information from fem to integration points using element based mapper
    #     mapper_parameters = KratosMultiphysics.Parameters("""{
    #         "mapper_type" : "nearest_element",
    #         "search_radius" : 1.0
    #     }""")
    #     mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, mapper_parameters )
    #     mapper.Map( KratosShape.SHAPE_CHANGE, KratosShape.SHAPE_CHANGE )


    #     # Assign shape updates to knot intersections
    #     list_of_nodes = [node for node in destination_mdpa.Nodes]
    #     list_of_updates = [node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE) for node in destination_mdpa.Nodes]

    #     for face in self.cad_model.GetByType('BrepFace'):
    #         face_key = face.Key()
    #         u_list = list(set(relevant_parameters[face_key]["u"]))
    #         v_list = list(set(relevant_parameters[face_key]["v"]))
    #         u_list.sort()
    #         v_list.sort()

    #         print(u_list)
    #         print(v_list)

    #         for v_itr in range(len(v_list)-2):

    #             print(v_itr)

    #             v0 = v_list[v_itr]
    #             v1 = v_list[v_itr+1]
    #             v2 = v_list[v_itr+2]

    #             for u_itr in range(len(u_list)-2):

    #                 u0 = u_list[u_itr]
    #                 u1 = u_list[u_itr+1]
    #                 u2 = u_list[u_itr+2]

    #                 print((face_key, u0, v0))

    #                 update00 = list_of_updates[ id_to_itr[(face_key, u0, v0)] ]
    #                 update10 = list_of_updates[ id_to_itr[(face_key, u1, v0)] ]
    #                 update20 = list_of_updates[ id_to_itr[(face_key, u2, v0)] ]
    #                 update01 = list_of_updates[ id_to_itr[(face_key, u0, v1)] ]
    #                 update02 = list_of_updates[ id_to_itr[(face_key, u0, v2)] ]

    #                 # Compute DUpdateDu


    #                 # Compute DUpdateDu




    #                     # mapped_shape_change[face.Key()][(u,v)] = list_of_nodes[node_itr].GetSolutionStepValue(KratosShape.SHAPE_CHANGE)

    #         # for face, (u,v) in zip(list_of_faces, list_of_parameters):
    #         #     geometry = face.Data().Geometry()

    #         #     distance = geometry_data.PointAt(u, v) - node_target_position

    #         #     node.SetSolutionStepValue(KratosShape.FITTING_ERROR, distance.tolist())

    #         #     if np.linalg.norm(distance) > refinement_tolerance:

    #         #         start_span = an.Knots.UpperSpan(geometry_data.DegreeU(), geometry_data.KnotsU(), u)
    #         #         end_span = start_span+1
    #         #         u_interval = (geometry_data.KnotsU()[start_span], geometry_data.KnotsU()[end_span])
    #         #         intervals_along_u_to_refine[geometry.Key()].append(u_interval)

    #         #         start_span = an.Knots.UpperSpan(geometry_data.DegreeV(), geometry_data.KnotsV(), v)
    #         #         end_span = start_span+1
    #         #         v_interval = (geometry_data.KnotsV()[start_span], geometry_data.KnotsV()[end_span])
    #         #         intervals_along_v_to_refine[geometry.Key()].append(v_interval)


    #     self.__SaveCadModel("cad_after_prior_refinement.iga")

    #     print("> Prior refinement finished in" ,round( time.time()-start_time, 3 ), " s.")


    #     err
    #     return