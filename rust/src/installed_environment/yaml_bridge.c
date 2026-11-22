/* libyaml owns YAML syntax and emission. This bridge exposes events without
 * copying its private parser/emitter layouts into the Rust implementation. */
#include <stdlib.h>
#include <string.h>
#include <yaml.h>

typedef struct {
    int kind;
    const unsigned char *anchor, *tag, *value;
    size_t length;
    int style, implicit, quoted;
    size_t index, line, column;
} priva_yaml_event;

typedef struct {
    int kind;
    const char *context, *problem;
    size_t context_index, context_line, context_column;
    size_t problem_index, problem_line, problem_column;
    size_t offset;
    int value;
} priva_yaml_error;

typedef int (*event_callback)(void *, const priva_yaml_event *);
typedef void (*error_callback)(void *, const priva_yaml_error *);

int priva_yaml_parse(const unsigned char *input, size_t length,
                     event_callback callback, error_callback failure, void *context) {
    yaml_parser_t parser;
    if (!yaml_parser_initialize(&parser)) return -1;
    yaml_parser_set_input_string(&parser, input, length);
    int status = 1;
    for (;;) {
        yaml_event_t event;
        if (!yaml_parser_parse(&parser, &event)) {
            priva_yaml_error error = {
                parser.error, parser.context, parser.problem,
                parser.context_mark.index, parser.context_mark.line, parser.context_mark.column,
                parser.problem_mark.index, parser.problem_mark.line, parser.problem_mark.column,
                parser.problem_offset, parser.problem_value
            };
            failure(context, &error);
            status = 0;
            break;
        }
        priva_yaml_event value = {0};
        value.kind = event.type;
        value.index = event.start_mark.index;
        value.line = event.start_mark.line;
        value.column = event.start_mark.column;
        switch (event.type) {
            case YAML_ALIAS_EVENT: value.anchor = event.data.alias.anchor; break;
            case YAML_SCALAR_EVENT:
                value.anchor = event.data.scalar.anchor;
                value.tag = event.data.scalar.tag;
                value.value = event.data.scalar.value;
                value.length = event.data.scalar.length;
                value.style = event.data.scalar.style;
                value.implicit = event.data.scalar.plain_implicit;
                value.quoted = event.data.scalar.quoted_implicit;
                break;
            case YAML_SEQUENCE_START_EVENT:
                value.anchor = event.data.sequence_start.anchor;
                value.tag = event.data.sequence_start.tag;
                value.implicit = event.data.sequence_start.implicit;
                value.style = event.data.sequence_start.style;
                break;
            case YAML_MAPPING_START_EVENT:
                value.anchor = event.data.mapping_start.anchor;
                value.tag = event.data.mapping_start.tag;
                value.implicit = event.data.mapping_start.implicit;
                value.style = event.data.mapping_start.style;
                break;
            default: break;
        }
        int continued = callback(context, &value);
        int ended = event.type == YAML_STREAM_END_EVENT;
        yaml_event_delete(&event);
        if (!continued || ended) break;
    }
    yaml_parser_delete(&parser);
    return status;
}

void *priva_yaml_emitter_new(yaml_write_handler_t callback, void *context) {
    yaml_emitter_t *emitter = malloc(sizeof(*emitter));
    if (!emitter) return NULL;
    if (!yaml_emitter_initialize(emitter)) { free(emitter); return NULL; }
    yaml_emitter_set_output(emitter, callback, context);
    yaml_emitter_set_unicode(emitter, 0);
    yaml_emitter_set_indent(emitter, 2);
    yaml_emitter_set_width(emitter, 80);
    yaml_emitter_set_break(emitter, YAML_LN_BREAK);
    return emitter;
}

int priva_yaml_emit(void *handle, const priva_yaml_event *value) {
    yaml_emitter_t *emitter = handle;
    yaml_event_t event;
    int ready = 0;
    switch (value->kind) {
        case YAML_STREAM_START_EVENT: ready = yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING); break;
        case YAML_STREAM_END_EVENT: ready = yaml_stream_end_event_initialize(&event); break;
        case YAML_DOCUMENT_START_EVENT: ready = yaml_document_start_event_initialize(&event, NULL, NULL, NULL, 1); break;
        case YAML_DOCUMENT_END_EVENT: ready = yaml_document_end_event_initialize(&event, 1); break;
        case YAML_ALIAS_EVENT: ready = yaml_alias_event_initialize(&event, value->anchor); break;
        case YAML_SCALAR_EVENT:
            ready = yaml_scalar_event_initialize(&event, value->anchor, value->tag,
                value->value, (int)value->length, value->implicit, value->quoted, YAML_ANY_SCALAR_STYLE);
            break;
        case YAML_SEQUENCE_START_EVENT:
            ready = yaml_sequence_start_event_initialize(&event, value->anchor,
                (const unsigned char *)YAML_SEQ_TAG, 1, YAML_BLOCK_SEQUENCE_STYLE);
            break;
        case YAML_SEQUENCE_END_EVENT: ready = yaml_sequence_end_event_initialize(&event); break;
        case YAML_MAPPING_START_EVENT:
            ready = yaml_mapping_start_event_initialize(&event, value->anchor,
                (const unsigned char *)YAML_MAP_TAG, 1, YAML_BLOCK_MAPPING_STYLE);
            break;
        case YAML_MAPPING_END_EVENT: ready = yaml_mapping_end_event_initialize(&event); break;
        default: return 0;
    }
    if (!ready) return 0;
    return yaml_emitter_emit(emitter, &event);
}

void priva_yaml_emitter_delete(void *handle) {
    yaml_emitter_t *emitter = handle;
    yaml_emitter_delete(emitter);
    free(emitter);
}
