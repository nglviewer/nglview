import { Dict } from "@jupyter-widgets/base"
import { StructureRef } from "molstar/lib/mol-plugin-state/manager/structure/hierarchy-state"
import { PluginUIContext } from "molstar/lib/mol-plugin-ui/context"

export
function addRepresentation(plugin: PluginUIContext, params, strucIndex){
    var st: StructureRef = plugin.managers.structure.hierarchy.current.structures[strucIndex]
    console.log("Calling from addRepresentation", st, params, strucIndex)
    var components = st.components
    plugin.dataTransaction(async () => {
        for (const component of components) {
            await plugin.builders.structure.representation.addRepresentation(
                component.cell, params)
        }
     })
}

export
function removeRepresentation(){
}

export
function clearRepresentation(){
}
