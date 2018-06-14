namespace eval Geo::examples {

}

proc Geo::examples::Init { } {
    uplevel #0 [list source [file join $::Geo::dir examples TrussCantilever.tcl]]
}

proc Geo::examples::UpdateMenus { } {
    GiDMenu::InsertOption "Kratos" [list "---"] 8 PRE "" "" "" insertafter =
    GiDMenu::InsertOption "Kratos" [list "Truss cantilever" ] 8 PRE [list ::Geo::examples::TrussCantilever] "" "" insertafter =
    GiDMenu::UpdateMenus
}

Geo::examples::Init
