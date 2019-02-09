using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Serialization;
using System.IO;


namespace Quantization_Tool
{
    class Data_Output
    {
        private XmlSerializer State_XmlSerializer;
        private XmlSerializer Simulation_XmlSerializer;
        private XmlSerializer Output_XmlSerializer;

        public Data_Output()
        {
            State_XmlSerializer = new XmlSerializer(typeof(State_Variables));
            Simulation_XmlSerializer = new XmlSerializer(typeof(Simulation));
            Output_XmlSerializer = new XmlSerializer(typeof(Output_Variables));
        }

        // Writing to File

        public void Write_XML(string file_name, State_Variables state)
        {
            TextWriter textWriter = new StreamWriter(file_name);

            State_XmlSerializer.Serialize(textWriter, state);

            textWriter.Close();
        }

        public void Write_XML(string file_name, Simulation simulation)
        {
            TextWriter textWriter = new StreamWriter(file_name);

            Simulation_XmlSerializer.Serialize(textWriter, simulation);

            textWriter.Close();
        }

        public void Write_XML(string file_name, Output_Variables output)
        {
            TextWriter textWriter = new StreamWriter(file_name);

            Output_XmlSerializer.Serialize(textWriter, output);

            textWriter.Close();
        }

        // Reading from File

        public State_Variables Sate_XML(string file_name)
        {
            FileStream fileStream  = new FileStream(file_name, FileMode.Open);

            State_Variables state = (State_Variables)State_XmlSerializer.Deserialize(fileStream);

            return state;
        }

        public Simulation Simulation_XML(string file_name)
        {
            FileStream fileStream = new FileStream(file_name, FileMode.Open);

            Simulation sim = (Simulation)Simulation_XmlSerializer.Deserialize(fileStream);

            return sim;
        }

        public Output_Variables Output_XML(string file_name)
        {
            FileStream fileStream = new FileStream(file_name, FileMode.Open);

            Output_Variables output = (Output_Variables)Output_XmlSerializer.Deserialize(fileStream);

            return output;
        }

    }
}
