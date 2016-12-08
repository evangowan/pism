// Copyright (C) 2011, 2013, 2014, 2015, 2016 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PCFACTORY_H_
#define _PCFACTORY_H_

#include <memory>
#include <map>

#include "base/util/IceGrid.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_options.hh"

namespace pism {

class Config;

template <class Model, class Modifier>
class PCFactory {
public:

  // virtual base class that allows storing different model creators
  // in the same dictionary
  class ModelCreator {
  public:
    virtual Model* create(IceGrid::ConstPtr g) = 0;
    virtual ~ModelCreator() {}
  };

  // Creator for a specific model class M.
  template <class M>
  class SpecificModelCreator : public ModelCreator {
  public:
    M* create(IceGrid::ConstPtr g) {
      return new M(g);
    }
  };

  // virtual base class that allows storing different modifier
  // creators in the same dictionary
  class ModifierCreator {
  public:
    virtual Modifier* create(IceGrid::ConstPtr g, Model* input) = 0;
    virtual ~ModifierCreator() {}
  };

  // Creator for a specific modifier class M.
  template <class M>
  class SpecificModifierCreator : public ModifierCreator {
  public:
    M* create(IceGrid::ConstPtr g, Model* input) {
      return new M(g, input);
    }
  };

  typedef std::shared_ptr<ModelCreator> ModelCreatorPtr;
  typedef std::shared_ptr<ModifierCreator> ModifierCreatorPtr;

  PCFactory<Model,Modifier>(IceGrid::ConstPtr g)
  : m_grid(g) {}
  ~PCFactory<Model,Modifier>() {}

  //! Sets the default type name.
  void set_default(std::string name) {
    if (m_models.find(name) == m_models.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "type %s is not registered", name.c_str());
    } else {
      m_default_type = name;
    }
  }

  //! Creates a boundary model. Processes command-line options.
  Model* create() {
    std::string model_list, modifier_list, description;
    Model* result = NULL;

    // build a list of available models:
    auto k = m_models.begin();
    model_list = "[" + (k++)->first;
    for (; k != m_models.end(); k++) {
      model_list += ", " + k->first;
    }
    model_list += "]";

    // build a list of available modifiers:
    auto p = m_modifiers.begin();
    modifier_list = "[" + (p++)->first;
    for (; p != m_modifiers.end(); p++) {
      modifier_list += ", " + p->first;
    }
    modifier_list += "]";

    description =  "Sets up the PISM " + m_option + " model. Available models: " + model_list +
      " Available modifiers: " + modifier_list;

    // Get the command-line option:
    options::StringList choices("-" + m_option, description, m_default_type);

    // the first element has to be an *actual* model (not a modifier), so we
    // create it:
    auto j = choices->begin();

    if (m_models.find(*j) == m_models.end()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s model \"%s\" is not available.\n"
                                    "Available models:    %s\n"
                                    "Available modifiers: %s",
                                    m_option.c_str(), j->c_str(),
                                    model_list.c_str(), modifier_list.c_str());
    }

    result = m_models[*j]->create(m_grid);

    ++j;

    // process remaining arguments:
    while (j != choices->end()) {
      if (m_modifiers.find(*j) == m_modifiers.end()) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s modifier \"%s\" is not available.\n"
                                      "Available modifiers: %s",
                                      m_option.c_str(), j->c_str(), modifier_list.c_str());
      }

      result =  m_modifiers[*j]->create(m_grid, result);

      ++j;
    }

    return result;
  }

  //! Adds a boundary model to the dictionary.
  template <class M>
  void add_model(const std::string &name) {
    m_models[name] = ModelCreatorPtr(new SpecificModelCreator<M>);
  }

  template <class M>
  void add_modifier(const std::string &name) {
    m_modifiers[name] = ModifierCreatorPtr(new SpecificModifierCreator<M>);
  }

  //! Removes a boundary model from the dictionary.
  void remove_model(const std::string &name) {
    m_models.erase(name);
  }

  void remove_modifier(const std::string &name) {
    m_modifiers.erase(name);
  }

  //! Clears the dictionary.
  void clear_models() {
    m_models.clear();
  }

  void clear_modifiers() {
    m_modifiers.clear();
  }
protected:
  std::string m_default_type, m_option;
  std::map<std::string, ModelCreatorPtr> m_models;
  std::map<std::string, ModifierCreatorPtr> m_modifiers;
  IceGrid::ConstPtr m_grid;
};

} // end of namespace pism

#endif /* _PCFACTORY_H_ */
