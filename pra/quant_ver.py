import quant_class as qc
# Enter your events here as follows
# Basic failure events of a component as such
# "<id>" : qc.event("<description>",<failure probability for on component>,"<consequence from list>",count=<number of components>,redundancy=<number of redudant components>)
# Failure event that require a combination of events as
# "<id>" : qc.comp("<description>","<consequence from list>",<list of id's contributing to event>)
event_dict ={
    "e1" : qc.event("test1",0.1,"LOC",count= 3),
    "e2" : qc.event("test1",0.1,"LOC",count= 3),
    "e3" : qc.event("test1",0.1,"LOM",count= 3,redundancy = 1)
}
compined_dict = {
    "c1" : qc.comp("test1","LOV",["e3","e1","e2"]),
    "c2" : qc.comp("test1","LOM",["c1","e1"])

}
test = qc.PRA(event_dict,compined_dict,["LOM","LOV","LOC"])
print(test.proability)