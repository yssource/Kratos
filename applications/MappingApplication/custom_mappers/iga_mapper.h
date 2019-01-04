//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_IGA_MAPPER_H_INCLUDED )
#define  KRATOS_IGA_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"
#include "custom_searching/interface_communicator.h"
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_flags.h"
#include "custom_utilities/mapper_local_system.h"


namespace Kratos
{
///@name Kratos Classes
///@{


class IgaLocalSystem : public MapperLocalSystem
{
public:

    explicit IgaLocalSystem(Condition* pIgaMapperCondition,
        NodePointerType pNode) : mpIgaMapperCondition(pIgaMapperCondition), mpNode(pNode) {}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_ERROR << "For now this function is not supposed to be used!" << std::endl;
    }

    /// Turn back information as a string.
    std::string PairingInfo(const int EchoLevel, const int CommRank) const override;

private:
    Condition* mpIgaMapperCondition;
    NodePointerType mpNode;
};

template<class TSparseSpace, class TDenseSpace>
class IgaMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IgaMapper
    KRATOS_CLASS_POINTER_DEFINITION(IgaMapper);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef Kratos::unique_ptr<InterfaceCommunicator> InterfaceCommunicatorPointerType;
    typedef typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef Kratos::unique_ptr<TMappingMatrixType> TMappingMatrixUniquePointerType;

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    IgaMapper(ModelPart& rModelPartIga,
              ModelPart& rModelPartFem)
                         : mrModelPartIga(rModelPartIga),
                           mrModelPartFem(rModelPartFem) {}

    IgaMapper(ModelPart& rModelPartIga,
              ModelPart& rModelPartFem,
              Parameters JsonParameters)
                : mrModelPartIga(rModelPartIga),
                  mrModelPartFem(rModelPartFem),
                  mMapperSettings(JsonParameters)
    {
        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartIga);
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartFem);

        ValidateInput(mMapperSettings);
        InitializeInterfaceCommunicator();
    }

    /// Destructor.
    ~IgaMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartIga,
                                  ModelPart& rModelPartFem,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<IgaMapper<TSparseSpace, TDenseSpace>>(
            rModelPartIga,
            rModelPartFem,
            JsonParameters);
    }

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        // // Set the Flags according to the type of remeshing
        // if (MappingOptions.Is(MapperFlags::REMESHED)) {
        //     InitializeInterface(MappingOptions);
        // }
        // else {
        //     BuildMappingMatrix(MappingOptions);
        // }

        // if (mpInverseMapper) {
        //     mpInverseMapper->UpdateInterface(MappingOptions,
        //                                      SearchRadius);
        // }
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        // if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
        //     GetInverseMapper()->Map(rOriginVariable, rDestinationVariable, MappingOptions);
        // }
        // else {
        //     MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        // }
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        // if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
        //     GetInverseMapper()->Map(rOriginVariable, rDestinationVariable, MappingOptions);
        // }
        // else {
        //     MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        // }
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        // if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
        //     MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        // }
        // else {
        //     GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        // }
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        // if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
        //     MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        // }
        // else {
        //     GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        // }
    }

    ///@}
    ///@name Access
    ///@{

    TMappingMatrixType* pGetMappingMatrix() override
    {
        return mpMappingMatrix.get();
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IgaMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }
protected:

   /**
    * @brief Initializing the Mapper
    * This has to be called in the constructor of the
    * derived classes, since it involves calls to
    * pure virtual functions
    */
    void Initialize()
    {
        InitializeInterface();
    }

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPartIga;
    ModelPart& mrModelPartFem;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    TMappingMatrixUniquePointerType mpMappingMatrix;

    MapperLocalSystemPointerVector mMapperLocalSystems;

    InterfaceCommunicatorPointerType mpIntefaceCommunicator;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;

    ///@}
    ///@name Private Operations
    ///@{

    void ValidateInput(Parameters AllMapperSettings);

    void ValidateParameters(Parameters AllMapperSettings)
    {
        Parameters default_settings = Parameters( R"({
            "search_radius"            : -1.0,
            "search_iterations"        : 3,
            "echo_level"               : 0
        })");

        AllMapperSettings.ValidateAndAssignDefaults(default_settings);
    }

    void InitializeInterfaceCommunicator();

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void BuildMappingMatrix(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartIga.GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartFem.GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<double>& rOriginVariable,
                              const Variable<double>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
                     const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<array_1d<double, 3>>& rOriginVariable,
                              const Variable<array_1d<double, 3>>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void PrintPairingInfo(const int EchoLevel);

    // // functions for customizing the behavior of this Mapper
    // virtual void CreateMapperLocalSystems(
    //     const Communicator& rModelPartCommunicator,
    //     std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) = 0;

    // virtual MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const = 0;

    ///@}
    ///@name Private  Access
    ///@{

    MapperUniquePointerType& GetInverseMapper()
    {
        if (!mpInverseMapper) {
            InitializeInverseMapper();
        }
        return mpInverseMapper;
    }

    void InitializeInverseMapper()
    {
        mpInverseMapper = this->Clone(mrModelPartFem,
                                      mrModelPartIga,
                                      mMapperSettings);
    }

}; // Class IgaMapper

}  // namespace Kratos.

#endif // KRATOS_IGA_MAPPER_H_INCLUDED  defined
