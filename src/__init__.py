from chimerax.core.toolshed import BundleAPI

class _MyAPI(BundleAPI):
    api_version = 1
    
    @staticmethod
    def start_tool(session, bi, ti):
        from .tool import PDBModifierTool
        return PDBModifierTool(session, ti.name)

bundle_api = _MyAPI()
