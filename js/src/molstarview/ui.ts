import { PluginLayoutControlsDisplay } from 'molstar/lib/mol-plugin/layout'

export
const basicSpec = {
            behaviors: [],
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: true,
                    controlsDisplay: 'reactive' as PluginLayoutControlsDisplay,
                    regionState: {
                        // bottom: 'full',
                        //left: 'collapsed',
                        //right: 'hidden',
                        // top: 'full',
                    }
                },
            },
}
