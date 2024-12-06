import { PluginUIContext } from "molstar/lib/mol-plugin-ui/context"

export
function addRepresentation(plugin: PluginUIContext, params, modelIndex){
    var st = plugin.managers.structure.hierarchy.current.structures[modelIndex]
    console.log("Calling from addRepresentation", st, params, modelIndex)
    var components = st.components
    plugin.dataTransaction(async () => {
        for (const component of components) {
            await plugin.builders.structure.representation.addRepresentation(
                component.cell, params)
        }
     })
}

export
function removeRepresentation(plugin: PluginUIContext, modelIndex){
}

export
function clearRepresentation(){
}
