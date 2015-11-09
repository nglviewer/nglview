/**
 * @file  MDSrv
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */


var MdsrvDatasource = function( baseUrl ){

    baseUrl = baseUrl || "";

    this.getListing = function( path ){
        path = path || "";
        var url = baseUrl + "dir/" + path;
        return NGL.autoLoad( url, {
            ext: "json", noWorker: true
        } ).then( function( jsonData ){
            return {
                path: path,
                data: jsonData.data
            };
        } );
    };

    this.getUrl = function( src ){
        var info = NGL.getFileInfo( src );
        return baseUrl + "file/" + info.path;
    };

    this.getNumframesUrl = function( src ){
        var info = NGL.getFileInfo( src );
        return baseUrl + "traj/numframes/" + info.path;
    };

    this.getFrameUrl = function( src, frameIndex ){
        var info = NGL.getFileInfo( src );
        return baseUrl + "traj/frame/" + frameIndex + "/" + info.path;
    };

    this.getFrameParams = function( src, atomIndices ){
        var info = NGL.getFileInfo( src );
        return "atomIndices=" + atomIndices.join(";");
    };

    this.getPathUrl = function( src, atomIndex ){
        var info = NGL.getFileInfo( src );
        return baseUrl + "traj/path/" + atomIndex + "/" + info.path;
    };

};
