
export
function addRepresentation(plugin, params, modelIndex){
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
function removeRepresentation(){
}

export 
function clearRepresentation(){
}
