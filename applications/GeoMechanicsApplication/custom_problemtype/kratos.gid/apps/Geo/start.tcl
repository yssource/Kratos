namespace eval ::Geo {
    # Variable declaration
    variable dir
    variable attributes
    variable kratos_name
}

proc ::Geo::Init { } {
    # Variable initialization
    variable dir
    variable attributes
    variable kratos_name
    
    set dir [apps::getMyDir "Geo"]
    set attributes [dict create]
    
    # Allow to open the tree
    set ::spdAux::TreeVisibility 1
    
    # Intervals 
    dict set attributes UseIntervals 1
    if {$::Kratos::kratos_private(DevMode) eq "dev"} {dict set attributes UseIntervals 1}
    
    set kratos_name GeoMechanicsApplication
    
    LoadMyFiles
}

proc ::Geo::LoadMyFiles { } {
    variable dir
    
    uplevel #0 [list source [file join $dir examples examples.tcl]]
    uplevel #0 [list source [file join $dir xml GetFromXML.tcl]]
    uplevel #0 [list source [file join $dir write write.tcl]]
    uplevel #0 [list source [file join $dir write writeProjectParameters.tcl]]
    uplevel #0 [list source [file join $dir postprocess formfinding.tcl]]
}

proc ::Geo::CustomToolbarItems { } {
    Kratos::ToolbarAddItem "Example" "example.png" [list -np- ::Geo::examples::TrussCantilever] [= "Example\nTruss cantilever"]   
}

proc ::Geo::CustomMenus { } {
    Geo::examples::UpdateMenus

    GiDMenu::InsertOption "Kratos" [list "---"] 8 PRE "" "" "" insertafter =
    GiDMenu::InsertOption "Kratos" [list "Formfinding - Update geometry" ] end POST [list ::Geo::Formfinding::UpdateGeometry] "" "" insert =
    GiDMenu::UpdateMenus
}

proc ::Geo::GetAttribute {name} {
    variable attributes
    set value ""
    if {[dict exists $attributes $name]} {set value [dict get $attributes $name]}
    return $value
}

proc ::Geo::BeforeMeshGeneration { size } { 
    foreach group [GiD_Groups list] {
        GiD_AssignData condition relation_line_geo_mesh Lines {0} [GiD_EntitiesGroups get $group lines]
    }
}

::Geo::Init
