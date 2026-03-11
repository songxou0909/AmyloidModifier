from chimerax.core.toolshed import BundleAPI

class _MyAPI(BundleAPI):
    # This explicitly tells ChimeraX to use the modern argument structure
    api_version = 1
    
    @staticmethod
    def start_tool(session, bi, ti):
        # Imports the tool class from your tool.py file and launches it
        from .tool import PDBModifierTool
        return PDBModifierTool(session, ti.name)

bundle_api = _MyAPI()