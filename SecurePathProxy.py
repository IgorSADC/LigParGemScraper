class SecurePathProxy:
    """
    This class is a good way to change directories without having to remember to go back later.
    Use with 'with'!
    There's also a very nice method to apply any callable entity in all children folders recursivilly
    """
    def __init__(self, os_interfaced_object, path):
        self.os_interfaced_object = os_interfaced_object
        self.path = path
        self.in_folder_event = lambda *args: None
        self.out_folder_event = lambda *args: None
        
    def __enter__(self):
        self.start_path = self.os_interfaced_object.getcwd()
        #The reason that is on an excpetion is because the sftp client doesnt have the path library. So I can't check if the path exists on there
        try:
            self.os_interfaced_object.chdir(self.path)
        except:
            self.os_interfaced_object.mkdir(self.path)
            self.os_interfaced_object.chdir(self.path)
        self.on_changing_folder()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.os_interfaced_object.chdir(self.start_path)
        
    def on_changing_folder(self):
        """This is a nice callback to whatever you want. It's called everytime there's a folder change"""
        if not self.in_folder_event: return
        self.in_folder_event(self.path)
    
    def on_backing_to_root_folder(self):
        if not self.out_folder_event: return
        self.out_folder_event()
        
    #THAT METHOD FOR FINDING FOLDERS IS NOT TRUSTABLE. 
    #THE ONLY REASON I'M USING IT IS BECAUSE THE SFTP CLIENT DOESN'T HAVE A DEFINITION FOR isdir()    
    def apply_callable_entity_recursivily(self, callable_entity, *args, **kwargs):
        if not callable(callable_entity): raise ValueError("Please pass a callable entity. If you are using a class, try overriding the durder method __call__")
        dir_list = []
        if hasattr(self.os_interfaced_object, 'path'):
            dir_list =  [a for a in self.os_interfaced_object.listdir() if self.os_interfaced_object.path.isdir(a)]
        else:
            dir_list = [a for a in self.os_interfaced_object.listdir() if len(a.split('.')) ==1]
          
        callable_entity(current_path = self.path, os_interfaced_object= self.os_interfaced_object,*args, **kwargs)
        
        for d in dir_list:
            with SecurePathProxy(self.os_interfaced_object, d) as new_path:
                new_path.apply_callable_entity_recursivily(callable_entity, *args, **kwargs)